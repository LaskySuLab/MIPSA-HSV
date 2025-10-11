#!/usr/bin/env Rscript
# Tripartite network plots: Virus → Auto-Ab → Disease
# Reproduces HSV-1/HSV-2/CMV/EBV panels with your layout logic
suppressPackageStartupMessages({
  library(argparse); library(data.table); library(dplyr); library(tidyr)
  library(igraph); library(ggplot2); library(ggrepel); library(scales)
})

# ----------------------
# CLI
# ----------------------
parser <- ArgumentParser()
parser$add_argument("--human_dx", required=TRUE, help="rep_bin_fchange_FL_dx_com_glm_ann.txt")
parser$add_argument("--proid_dx", required=TRUE, help="rep_bin_fchange_HSV_dx_com_glm_ann.txt")
parser$add_argument("--human_qtl", required=TRUE, help="rep_hsv.sig (or equivalent QTL table)")
parser$add_argument("--species", nargs="+", required=TRUE,
                    help="One or more taxon_species, e.g. 'Human alphaherpesvirus 1' 'Human alphaherpesvirus 2'")
parser$add_argument("--p_thresh", type="double", default=0.05,
                    help="P-value threshold for edges (default 0.05)")
parser$add_argument("--outdir", required=TRUE, help="Output directory for PNGs")
parser$add_argument("--seed", type="integer", default=123, help="Seed for layout_with_fr")
args <- parser$parse_args()

dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)

# ----------------------
# Read + pre-clean
# ----------------------
human_dx <- fread(args$human_dx)
proid_dx <- fread(args$proid_dx)
human_qtl <- fread(args$human_qtl)

# Fallback gene symbols if missing
human_dx$gene_symbol  <- ifelse(is.na(human_dx$gene_symbol),  human_dx$UniProt_acc,  human_dx$gene_symbol)
proid_dx$gene_symbol  <- ifelse(is.na(proid_dx$gene_symbol),  proid_dx$UniProt_acc,  proid_dx$gene_symbol)
human_qtl$gene_symbol <- ifelse(is.na(human_qtl$gene_symbol), human_qtl$UniProt_acc, human_qtl$gene_symbol)

# Make disease labels pretty
human_dx$disease <- gsub("_", " ", human_dx$disease)
proid_dx$disease <- gsub("_", " ", proid_dx$disease)

# map long species → short alias for titles/filenames
alias_map <- c(
  "Human alphaherpesvirus 1" = "HSV1",
  "Human alphaherpesvirus 2" = "HSV2",
  "Human alphaherpesvirus 3" = "VZV",
  "Human betaherpesvirus 5"  = "CMV",
  "Human betaherpesvirus 6A" = "HHV6A",
  "Human betaherpesvirus 6B" = "HHV6B",
  "Human betaherpesvirus 7"  = "HHV7",
  "Human gammaherpesvirus 4" = "EBV",
  "Human gammaherpesvirus 8" = "HHV8"
)

# ----------------------
# Helper: build one network for a species
# ----------------------
plot_one_species <- function(species_name, seed=args$seed, p_thr=args$p_thresh) {

  # Filter QTLs to species
  qtl_sp <- human_qtl %>% filter(taxon_species == species_name)

  # Subset association tables to this species' markers
  # human_dx: var_id = antibody (Hu)
  # proid_dx: var_id = viral UniProt (Virus)
  human_dx_sub <- human_dx %>% filter(var_id %in% qtl_sp$antibody)
  proid_dx_sub <- proid_dx %>% filter(var_id %in% qtl_sp$var_id)

  # Edges (p < thresh): widths by -log10(p)
  edges_vh <- qtl_sp %>%
    filter(!is.na(P), P < p_thr) %>%
    transmute(from = var_id, to = gene_symbol, logp = -log10(P), etype = "Virus→Auto-Ab")

  edges_hd <- human_dx_sub %>%
    filter(!is.na(P), P < p_thr) %>%
    transmute(from = gene_symbol, to = disease, logp = -log10(P), etype = "Auto-Ab→Disease")

  edges_vd <- proid_dx_sub %>%
    filter(!is.na(P), P < p_thr) %>%
    transmute(from = var_id, to = disease, logp = -log10(P), etype = "Virus→Disease")

  edges_all <- bind_rows(edges_vh, edges_hd, edges_vd)
  if (nrow(edges_all) == 0L) {
    message("No edges after P<", p_thr, " for ", species_name, " — skipping.")
    return(invisible(NULL))
  }

  # Nodes + types
  node_names <- unique(c(edges_all$from, edges_all$to))
  nodes <- tibble(name = node_names) %>%
    mutate(
      ntype = case_when(
        name %in% qtl_sp$var_id ~ "Virus",
        name %in% qtl_sp$gene_symbol ~ "Auto-Ab",
        name %in% unique(c(human_dx_sub$disease, proid_dx_sub$disease)) ~ "Disease",
        TRUE ~ "Other"
      )
    )

  # Size encodings
  ha_pep <- edges_vh %>% count(to,  name = "n_pep") %>% rename(name = to)
  ha_dis <- edges_hd %>% count(from, name = "n_dis") %>% rename(name = from)
  ha_size <- full_join(ha_pep, ha_dis, by = "name") %>%
    mutate(across(c(n_pep, n_dis), ~replace_na(., 0)), ha_total = n_pep + n_dis)

  dis_v  <- edges_vd %>% count(to,  name = "n_from_viral") %>% rename(name = to)
  dis_h  <- edges_hd %>% count(to,  name = "n_from_human") %>% rename(name = to)
  dis_size <- full_join(dis_v, dis_h, by = "name") %>%
    mutate(across(c(n_from_viral, n_from_human), ~replace_na(., 0)),
           dis_total = n_from_viral + n_from_human)

  diseases_with_human <- unique(edges_hd$to)

  nodes <- nodes %>%
    left_join(ha_size, by = "name") %>%
    left_join(dis_size, by = "name") %>%
    mutate(
      node_fill = case_when(
        ntype == "Auto-Ab" ~ "skyblue",
        ntype == "Virus"   ~ "green",
        ntype == "Disease" & name %in% diseases_with_human ~ "red",
        ntype == "Disease" ~ "pink",
        TRUE ~ "grey70"
      ),
      node_size = case_when(
        ntype == "Auto-Ab" ~ rescale(replace_na(ha_total, 0),  to = c(8, 21)),
        ntype == "Disease" ~ rescale(replace_na(dis_total, 0), to = c(4, 12)),
        ntype == "Virus"   ~ 2.5,
        TRUE ~ 2.5
      )
    )

  # Graph + initial FR layout
  g <- graph_from_data_frame(edges_all, vertices = nodes, directed = TRUE)
  w <- edges_all$logp; w[!is.finite(w) | w <= 0] <- 1e-3
  set.seed(seed)
  lay <- layout_with_fr(g, weights = w)
  layout_df <- as.data.frame(lay); colnames(layout_df) <- c("x","y"); layout_df$name <- V(g)$name

  nodes_plot <- nodes %>% left_join(layout_df, by = "name")

  # Put disease nodes on a ring (evenly spaced)
  cx <- mean(nodes_plot$x[nodes_plot$ntype != "Disease"])
  cy <- mean(nodes_plot$y[nodes_plot$ntype != "Disease"])

  r_non_dis <- max(sqrt((nodes_plot$x[nodes_plot$ntype!="Disease"]-cx)^2 +
                          (nodes_plot$y[nodes_plot$ntype!="Disease"]-cy)^2), na.rm = TRUE)
  r_outer <- r_non_dis * 1.2

  dis_df <- nodes_plot %>%
    filter(ntype == "Disease") %>%
    mutate(theta = atan2(y - cy, x - cx)) %>%
    arrange(theta) %>%
    mutate(theta_even = seq(-pi, pi, length.out = n()))

  dis_df <- dis_df %>%
    mutate(
      x_adj = cx + r_outer * cos(theta_even),
      y_adj = cy + r_outer * sin(theta_even)
    ) %>%
    select(name, x_adj, y_adj)

  nodes_plot <- nodes_plot %>%
    left_join(dis_df, by = "name") %>%
    mutate(
      x_adj = ifelse(is.na(x_adj), x, x_adj),
      y_adj = ifelse(is.na(y_adj), y, y_adj)
    )

  # Virus nodes on a middle ring aiming toward their targets
  nodes_pos <- nodes_plot %>%
    mutate(x_use = ifelse(is.na(x_adj), x, x_adj),
           y_use = ifelse(is.na(y_adj), y, y_adj))

  cx <- mean(nodes_pos$x_use[nodes_pos$ntype == "Auto-Ab"], na.rm = TRUE)
  cy <- mean(nodes_pos$y_use[nodes_pos$ntype == "Auto-Ab"], na.rm = TRUE)

  r_auto_max <- max(sqrt((nodes_pos$x_use[nodes_pos$ntype == "Auto-Ab"] - cx)^2 +
                           (nodes_pos$y_use[nodes_pos$ntype == "Auto-Ab"] - cy)^2), na.rm = TRUE)
  r_mid <- r_auto_max + 0.3 * (r_outer - r_auto_max)

  targets_pos <- nodes_pos %>% select(name, x_use, y_use)
  viral_names <- nodes_pos %>% filter(ntype == "Virus") %>% pull(name)

  viral_angles <- edges_all %>%
    filter(from %in% viral_names) %>%
    select(from, to) %>%
    left_join(targets_pos, by = c("to" = "name")) %>%
    mutate(theta = atan2(y_use - cy, x_use - cx)) %>%
    group_by(from) %>%
    summarise(theta = atan2(mean(sin(theta), na.rm = TRUE), mean(cos(theta), na.rm = TRUE)),
              .groups = "drop") %>%
    arrange(theta) %>%
    mutate(theta = theta + (row_number() - mean(row_number())) * 0.02)

  viral_middle <- viral_angles %>%
    transmute(name  = from,
              x_mid = cx + r_mid * cos(theta),
              y_mid = cy + r_mid * sin(theta))

  nodes_plot <- nodes_plot %>%
    left_join(viral_middle, by = "name") %>%
    mutate(
      x_adj = ifelse(ntype == "Virus" & !is.na(x_mid), x_mid, x_adj),
      y_adj = ifelse(ntype == "Virus" & !is.na(y_mid), y_mid, y_adj)
    ) %>%
    select(-x_mid, -y_mid)

  # Edge coordinates
  edges_plot <- edges_all %>%
    left_join(nodes_plot %>% select(name, x_adj, y_adj), by = c("from" = "name")) %>%
    rename(x_from = x_adj, y_from = y_adj) %>%
    left_join(nodes_plot %>% select(name, x_adj, y_adj), by = c("to" = "name")) %>%
    rename(x_to = x_adj, y_to = y_adj) %>%
    mutate(order_rank = case_when(
      etype == "Virus→Disease"   ~ 1,
      etype == "Virus→Auto-Ab"   ~ 2,
      etype == "Auto-Ab→Disease" ~ 3,
      TRUE ~ 4
    )) %>%
    arrange(order_rank)

  # Plot
  title_alias <- ifelse(species_name %in% names(alias_map), alias_map[[species_name]], species_name)
  p <- ggplot() +
    geom_segment(
      data = edges_plot,
      aes(x = x_from, y = y_from, xend = x_to, yend = y_to,
          color = etype, linewidth = logp),
      alpha = 0.50
    ) +
    scale_color_manual(values = c("Virus→Disease"   = "palegreen1",
                                  "Virus→Auto-Ab"   = "royalblue2",
                                  "Auto-Ab→Disease" = "red"),
                       name = "Edge Type") +
    guides(color = guide_legend(override.aes = list(linewidth = 2.8))) +
    scale_linewidth_continuous(range = c(1, 3), name = "-log10(p)") +
    geom_point(
      data = nodes_plot,
      aes(x = x_adj, y = y_adj, size = node_size, fill = node_fill),
      shape = 21, color = "grey20", stroke = 0.2
    ) +
    scale_fill_identity(name = "Node Type") +
    scale_size_identity() +
    geom_text(
      data = nodes_plot %>% filter(ntype == "Auto-Ab"),
      aes(x = x_adj, y = y_adj, label = name),
      size = 4, color = "black", fontface = "bold"
    ) +
    geom_text_repel(
      data = nodes_plot %>% filter(ntype == "Disease"),
      aes(x = x_adj, y = y_adj, label = name),
      size = 3.5, color = 'black'
    ) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.text  = element_text(size = 10),
          legend.title = element_text(size = 10, face = "bold")) +
    labs(title = paste0(title_alias, " Network (Virus → Auto-Ab → Disease)"))

  outfile <- file.path(args$outdir, paste0("network_plot_", title_alias, ".png"))
  ggsave(outfile, plot = p, width = 10, height = 10, dpi = 300)
  message("Saved: ", outfile)
}

# ----------------------
# Run for each requested species
# ----------------------
for (sp in args$species) {
  plot_one_species(sp)
}
