#!/usr/bin/env Rscript
# Tripartite network plots: Virus → Auto-Ab → Disease
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(dplyr); library(tidyr); library(igraph)
  library(ggplot2); library(ggrepel); library(scales); library(grid)
})

option_list <- list(
  make_option("--human_dx", type="character", help="rep_bin_fchange_FL_dx_com_glm_ann.txt"),
  make_option("--proid_dx", type="character", help="rep_bin_fchange_HSV_dx_com_glm_ann.txt"),
  make_option("--rep_sig",  type="character", help="rep_hsv.sig"),
  make_option("--virus",    type="character", help="Virus species to plot (e.g., 'Human alphaherpesvirus 1')"),
  make_option("--out_png",  type="character", default="results/network_plot.png", help="Output PNG filename")
)
opt <- parse_args(OptionParser(option_list=option_list))

human_dx <- fread(opt$human_dx)
proid_dx <- fread(opt$proid_dx)
human_qtl <- fread(opt$rep_sig)

human_dx <- human_dx %>%
  filter(!disease %in% c("Asthma","Depression","Headache","Migraine","Opioid_use_disorder")) %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol), UniProt_acc, gene_symbol),
         disease = gsub("_"," ", disease))

proid_dx <- proid_dx %>%
  filter(!disease %in% c("Asthma","Depression","Headache","Migraine","Opioid_use_disorder")) %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol), UniProt_acc, gene_symbol),
         disease = gsub("_"," ", disease))

human_qtl <- human_qtl %>%
  filter(taxon_species == opt$virus) %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol), UniProt_acc, gene_symbol))

human_dx_sub <- human_dx %>% filter(var_id %in% human_qtl$antibody)
proid_dx_sub <- proid_dx %>% filter(var_id %in% human_qtl$var_id)

edges_vh <- human_qtl %>%
  filter(!is.na(P), P < 0.05) %>%
  transmute(from = var_id, to = gene_symbol, logp = -log10(P), etype = "Virus→Auto-Ab")

edges_hd <- human_dx_sub %>%
  filter(!is.na(P), P < 0.05) %>%
  transmute(from = gene_symbol, to = disease, logp = -log10(P), etype = "Auto-Ab→Disease")

edges_vd <- proid_dx_sub %>%
  filter(!is.na(P), P < 0.05) %>%
  transmute(from = var_id, to = disease, logp = -log10(P), etype = "Virus→Disease")

edges_all <- bind_rows(edges_vh, edges_hd, edges_vd)
if (!nrow(edges_all)) stop("No significant edges (p < 0.05).")

node_names <- unique(c(edges_all$from, edges_all$to))
nodes <- tibble::tibble(name = node_names) %>%
  mutate(
    ntype = case_when(
      name %in% human_qtl$var_id ~ "Virus",
      name %in% human_qtl$gene_symbol ~ "Auto-Ab",
      name %in% unique(c(human_dx_sub$disease, proid_dx_sub$disease)) ~ "Disease",
      TRUE ~ "Other"
    )
  )

ha_pep <- edges_vh %>% count(to,  name = "n_pep") %>% rename(name = to)
ha_dis <- edges_hd %>% count(from, name = "n_dis") %>% rename(name = from)
ha_size <- full_join(ha_pep, ha_dis, by = "name") %>% mutate(across(everything(), ~replace_na(.,0)), ha_total = n_pep + n_dis)

dis_v  <- edges_vd %>% count(to,  name = "n_from_viral") %>% rename(name = to)
dis_h  <- edges_hd %>% count(to,  name = "n_from_human") %>% rename(name = to)
dis_size <- full_join(dis_v, dis_h, by = "name") %>% mutate(across(everything(), ~replace_na(.,0)), dis_total = n_from_viral + n_from_human)

diseases_with_human <- unique(edges_hd$to)

nodes <- nodes %>%
  left_join(ha_size, by="name") %>%
  left_join(dis_size, by="name") %>%
  mutate(
    node_fill = case_when(
      ntype == "Auto-Ab" ~ "skyblue",
      ntype == "Virus" ~ "green",
      ntype == "Disease" & name %in% diseases_with_human ~ "red",
      ntype == "Disease" ~ "pink",
      TRUE ~ "grey70"
    ),
    node_size = case_when(
      ntype == "Auto-Ab" ~ scales::rescale(replace_na(ha_total,0), to=c(8,21)),
      ntype == "Disease" ~ scales::rescale(replace_na(dis_total,0), to=c(4,12)),
      ntype == "Virus"   ~ 2.5,
      TRUE ~ 2.5
    )
  )

g <- graph_from_data_frame(edges_all, vertices = nodes, directed = TRUE)
w <- edges_all$logp; w[!is.finite(w) | w <= 0] <- 1e-3
set.seed(123)
lay <- layout_with_fr(g, weights = w)
layout_df <- as.data.frame(lay); colnames(layout_df) <- c("x","y"); layout_df$name <- V(g)$name

nodes_plot <- nodes %>% left_join(layout_df, by="name")

# Put diseases on outer ring, space evenly; center by auto-Abs
cx <- mean(nodes_plot$x[nodes_plot$ntype != "Disease"]); cy <- mean(nodes_plot$y[nodes_plot$ntype != "Disease"])
r_non <- max(sqrt((nodes_plot$x[nodes_plot$ntype!="Disease"]-cx)^2 + (nodes_plot$y[nodes_plot$ntype!="Disease"]-cy)^2), na.rm = TRUE)
r_outer <- r_non * 1.2
dis_df <- nodes_plot %>% filter(ntype=="Disease") %>%
  mutate(theta = atan2(y - cy, x - cx)) %>% arrange(theta) %>% mutate(theta_even = seq(-pi, pi, length.out = n())) %>%
  transmute(name, x_adj = cx + r_outer * cos(theta_even), y_adj = cy + r_outer * sin(theta_even))
nodes_plot <- nodes_plot %>% left_join(dis_df, by="name") %>%
  mutate(x_adj = ifelse(is.na(x_adj), x, x_adj), y_adj = ifelse(is.na(y_adj), y, y_adj))

# Edge coords
edges_plot <- edges_all %>% left_join(nodes_plot[, c("name","x_adj","y_adj")], by=c("from"="name")) %>% rename(x_from=x_adj, y_from=y_adj) %>%
  left_join(nodes_plot[, c("name","x_adj","y_adj")], by=c("to"="name")) %>% rename(x_to=x_adj, y_to=y_adj) %>%
  mutate(order_rank = dplyr::case_when(etype=="Virus→Disease" ~ 1, etype=="Virus→Auto-Ab" ~ 2, etype=="Auto-Ab→Disease" ~ 3, TRUE~4)) %>% arrange(order_rank)

p <- ggplot() +
  geom_segment(data=edges_plot, aes(x=x_from, y=y_from, xend=x_to, yend=y_to, color=etype, linewidth=logp), alpha=0.5) +
  scale_color_manual(values=c("Virus→Disease"="palegreen1","Virus→Auto-Ab"="royalblue2","Auto-Ab→Disease"="red"), name="Edge Type") +
  guides(color=guide_legend(override.aes = list(linewidth=2.8))) +
  scale_linewidth_continuous(range=c(1,3), name="-log10(p)") +
  geom_point(data=nodes_plot, aes(x=x_adj, y=y_adj, size=node_size, fill=node_fill), shape=21, color="grey20", stroke=0.2) +
  scale_fill_identity() + scale_size_identity() +
  geom_text(data=nodes_plot %>% filter(ntype=="Auto-Ab"), aes(x=x_adj, y=y_adj, label=name), size=4, fontface="bold") +
  geom_text_repel(data=nodes_plot %>% filter(ntype=="Disease"), aes(x=x_adj, y=y_adj, label=name), size=3.5) +
  theme_void() + theme(legend.position="bottom") +
  labs(title = paste0(opt$virus, " network (Virus → Auto-Ab → Disease)"))
ggsave(opt$out_png, p, width=10, height=10, dpi=300)
