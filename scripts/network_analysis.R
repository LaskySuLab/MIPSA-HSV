#!/usr/bin/env Rscript
# Tripartite network plots: Virus -> Auto-Ab -> Disease
# Layout: Concentric Layers (Center=Virus, Middle=Auto-Ag, Outer=Disease)

suppressPackageStartupMessages({
  library(optparse); library(data.table); library(dplyr); library(tidyr); library(igraph);
  library(ggplot2); library(ggrepel); library(scales); library(grid)
  source("R/utils_common.R"); source("R/utils_qtl.R")
})

option_list <- list(
  make_option("--human_dx_inc", type="character", help="Human Incidence file"),
  make_option("--human_dx_pre", type="character", help="Human Prevalence file"),
  make_option("--virus_dx_inc", type="character", help="Virus Incidence file"),
  make_option("--virus_dx_pre", type="character", help="Virus Prevalence file"),
  make_option("--sig_ab",       type="character", help="Significant Antibody file"),
  make_option("--out_dir",      type="character", default="results/Network", help="Output dir")
)

opt <- parse_args(OptionParser(option_list=option_list))
if(!dir.exists(opt$out_dir)) dir.create(opt$out_dir, recursive = TRUE)

# 1. Define Disease or LAB Categories
# ----------------------
dx_lookup <- stack(dx_categories)
names(dx_lookup) <- c("Disease", "Category")
dx_lookup$Disease <- gsub("_", " ", dx_lookup$Disease)

lab_lookup <- stack(lab_categories)
names(lab_lookup) <- c("LAB", "Category")
lab_lookup$LAB <- gsub("_most_recent", "", lab_lookup$LAB)
lab_lookup$LAB <- gsub("_", " ", lab_lookup$LAB)

# 2. Load Data
# ----------------------
cat("Loading data...\n")
human_dx.inc <- fread(opt$human_dx_inc)
virus_dx.inc <- fread(opt$virus_dx_inc)
human_dx.pre <- fread(opt$human_dx_pre)
virus_dx.pre <- fread(opt$virus_dx_pre)

human_lab <- fread(opt$human_lab)
virus_lab <- fread(opt$virus_lab)
sig_ab       <- fread(opt$sig_ab)

# Apply Label Mapping
virus_dx.inc$taxon_species <- species_label_map(virus_dx.inc$taxon_species)
virus_dx.pre$taxon_species <- species_label_map(virus_dx.pre$taxon_species)

virus_lab$taxon_species <- species_label_map(virus_lab$taxon_species)
sig_ab$taxon_species       <- species_label_map(sig_ab$taxon_species)

human_dx.inc$gene_symbol <- gsub(" \\s*\\([^\\)]+\\)", "",  human_dx.inc$gene_symbol)
human_dx.pre$gene_symbol <- gsub(" \\s*\\([^\\)]+\\)", "",  human_dx.pre$gene_symbol)
human_dx.inc$gene_symbol <- gsub("GPRASP3,ARMCX5-GPRASP2", "GPRASP3",  human_dx.inc$gene_symbol)
human_dx.pre$gene_symbol <- gsub("GPRASP3,ARMCX5-GPRASP2", "GPRASP3",  human_dx.pre$gene_symbol)

sig_ab$gene_symbol       <- gsub(" \\s*\\([^\\)]+\\)", "",  sig_ab$gene_symbol)
sig_ab$gene_symbol       <- gsub("GPRASP3,ARMCX5-GPRASP2", "GPRASP3",  sig_ab$gene_symbol)

human_lab$gene_symbol <- gsub(" \\s*\\([^\\)]+\\)", "",  human_lab$gene_symbol)


# 3. Core Plotting Function for diseases
# ----------------------
run_network_analysis <- function(target_virus, outcome_type, human_dx, virus_dx, sig_ab_full) {
  
  cat(sprintf("Processing: %s - %s\n", target_virus, outcome_type))
  
  # --- Step A: Filter Data (NOMINAL P < 0.05) ---
  
  # 1. Virus -> Auto-Ag
  edges_vh <- sig_ab_full %>%
    filter(taxon_species == target_virus) %>%
    filter(!is.na(P), P < 0.05) %>%
    transmute(from = var_id, to = gene_symbol, 
              pval = P, fdr = P.adj, 
              etype = "Virus→Auto-Ag")
  
  if(nrow(edges_vh) == 0) {
    cat(sprintf("  Skipping %s: No nominal Virus-Ab signals.\n", target_virus))
    return(NULL)
  }
  
  relevant_auto_ags <- unique(edges_vh$to)
  
  # 2. Auto-Ag -> Disease
  edges_hd <- human_dx %>%
    filter(gene_symbol %in% relevant_auto_ags) %>%
    filter(!is.na(P_all), P_all < 0.05) %>%
    transmute(from = gene_symbol, to = disease, 
              pval = P_all, fdr = p.adj_all, 
              etype = "Auto-Ag→Disease")
  
  # 3. Virus -> Disease
  edges_vd <- virus_dx %>%
    filter(taxon_species == target_virus) %>%
    filter(!is.na(P_all), P_all < 0.05) %>%
    transmute(from = UniProt_acc, to = disease, 
              pval = P_all, fdr = p.adj_all, 
              etype = "Virus→Disease")
  
  edges_all <- bind_rows(edges_vh, edges_hd, edges_vd)
  
  if(nrow(edges_all) == 0) return(NULL)
  
  # --- Step B: Edge Classification ---
  
  edges_all <- edges_all %>%
    mutate(
      logFDR = -log10(fdr),
      edge_class = case_when(
        etype == "Auto-Ag→Disease" & fdr < 0.05 ~ "Sig_AutoAg",
        (etype == "Virus→Auto-Ag" | etype == "Virus→Disease") & fdr < 0.05 ~ "Sig_Virus",
        TRUE ~ "Nominal"
      ),
      draw_order = case_when(
        edge_class == "Nominal" ~ 1,
        edge_class == "Sig_Virus" ~ 2,
        edge_class == "Sig_AutoAg" ~ 3
      )
    ) %>%
    arrange(draw_order)
  
  # --- Step C: Node definitions & Sizing ---
  
  node_names <- unique(c(edges_all$from, edges_all$to))
  nodes <- tibble::tibble(name = node_names) %>%
    mutate(
      ntype = case_when(
        name %in% sig_ab_full$var_id ~ "Virus", 
        name %in% sig_ab_full$gene_symbol ~ "Auto-Ag",
        name %in% unique(c(human_dx$disease, virus_dx$disease)) ~ "Disease",
        TRUE ~ "Other"
      )
    )
  
  # Sizing logic
  ha_pep <- edges_vh %>% count(to, name = "n_pep") %>% rename(name = to)
  ha_dis <- edges_hd %>% count(from, name = "n_dis") %>% rename(name = from)
  ha_size <- full_join(ha_pep, ha_dis, by = "name") %>% 
    mutate(across(c(n_pep, n_dis), ~replace_na(.,0)), ha_total = n_pep + n_dis)
  
  dis_v <- edges_vd %>% count(to, name = "n_from_viral") %>% rename(name = to)
  dis_h <- edges_hd %>% count(to, name = "n_from_human") %>% rename(name = to)
  dis_size <- full_join(dis_v, dis_h, by = "name") %>% 
    mutate(across(c(n_from_viral, n_from_human), ~replace_na(.,0)), dis_total = n_from_viral + n_from_human)
  
  edges_hd.sig = subset(edges_hd, edges_hd$fdr <0.05)
  diseases_with_human <- unique(edges_hd.sig$to)
  edges_vd.sig = subset(edges_vd, edges_vd$fdr <0.05)
  diseases_with_virus <- unique(edges_vd.sig$to)
  
  nodes_final <- nodes %>%
    left_join(ha_size, by="name") %>%
    left_join(dis_size, by="name") %>%
    mutate(
      node_fill = case_when(
        ntype == "Auto-Ag" ~ "skyblue",
        ntype == "Virus" ~ "green",
        ntype == "Disease" & name %in% c(diseases_with_human,diseases_with_virus) ~ "red",
        ntype == "Disease" ~ "pink",
        TRUE ~ "grey70"
      ),
      node_size = case_when(
        ntype == "Auto-Ag" ~ scales::rescale(replace_na(ha_total, 0), to=c(5, 15)), 
        ntype == "Disease" ~ scales::rescale(replace_na(dis_total, 0), to=c(4, 10)),
        ntype == "Virus"   ~ 2.5,
        TRUE ~ 2.5
      )
    )
  nodes_final$node_size[is.na(nodes_final$node_size)] <- 3
  
  # --- Step D: Concentric Layout Calculation ---
  # CONFIGURATION:
  # 1. Inner Core (R=2): Virus Peptides
  # 2. Middle Ring (R=6): Auto-Ags
  # 3. Outer Ring (R=10): Diseases
  
  r_inner  <- 2
  r_middle <- 6
  r_outer  <- 10
  
  nodes_plot <- nodes_final
  
  # D.1: Diseases (Outer Ring - Grouped)
  dis_subset <- nodes_plot %>% filter(ntype == "Disease")
  if(nrow(dis_subset) > 0) {
    dis_df <- dis_subset %>% left_join(dx_lookup, by = c("name" = "Disease")) %>%
      arrange(Category, name)
    
    dis_df$cat_id <- as.numeric(factor(dis_df$Category, levels=unique(dis_df$Category)))
    gap_size <- 2
    dis_df <- dis_df %>%
      mutate(
        row_idx = row_number(),
        total_slots = n() + (max(cat_id) * gap_size),
        pos_idx = row_idx + (cat_id * gap_size),
        theta = (pos_idx / max(pos_idx)) * 2 * pi,
        x_adj = r_outer * cos(theta),
        y_adj = r_outer * sin(theta)
      )
    
    nodes_plot <- nodes_plot %>% 
      left_join(dis_df %>% select(name, x_adj, y_adj), by="name")
  } 
  
  # D.2: Auto-Ags (Middle Ring - Radius 6) -> MOVED FROM CENTER TO MIDDLE
  ag_subset <- nodes_plot %>% filter(ntype == "Auto-Ag")
  if(nrow(ag_subset) > 0) {
    theta_ag <- seq(0, 2*pi, length.out = nrow(ag_subset) + 1)[1:nrow(ag_subset)]
    ag_df <- tibble(
      name = ag_subset$name,
      x_adj = r_middle * cos(theta_ag),
      y_adj = r_middle * sin(theta_ag)
    )
    
    nodes_plot <- nodes_plot %>% 
      rows_update(ag_df, by="name")
  }
  
  # D.3: Virus Peptides (Inner Core - Radius 2) -> MOVED FROM MIDDLE TO CENTER
  vir_subset <- nodes_plot %>% filter(ntype == "Virus")
  if(nrow(vir_subset) > 0) {
    theta_vir <- seq(0, 2*pi, length.out = nrow(vir_subset) + 1)[1:nrow(vir_subset)]
    vir_df <- tibble(
      name = vir_subset$name,
      x_adj = r_inner * cos(theta_vir),
      y_adj = r_inner * sin(theta_vir)
    )
    
    nodes_plot <- nodes_plot %>% 
      rows_update(vir_df, by="name")
  }
  
  # Fallback for unconnected nodes
  nodes_plot$x_adj[is.na(nodes_plot$x_adj)] <- 0
  nodes_plot$y_adj[is.na(nodes_plot$y_adj)] <- 0
  
  # Map coordinates to edges
  edges_plot <- edges_all %>%
    left_join(nodes_plot[, c("name","x_adj","y_adj")], by=c("from"="name")) %>% rename(x_from=x_adj, y_from=y_adj) %>%
    left_join(nodes_plot[, c("name","x_adj","y_adj")], by=c("to"="name")) %>% rename(x_to=x_adj, y_to=y_adj)
  
  # --- Step E: Plotting ---
  set.seed(2026)
  p <- ggplot() +
    # Edges
    geom_segment(data=edges_plot, 
                 aes(x=x_from, y=y_from, xend=x_to, yend=y_to, 
                     color=edge_class, linewidth=edge_class, alpha=edge_class), 
                 lineend = "round") +
    
    # Scales
    scale_color_manual(values=c("Sig_Virus"  = "darkgreen", 
                                "Sig_AutoAg" = "red", 
                                "Nominal"    = "grey50"), name="Significance") +
    
    scale_linewidth_manual(values=c("Sig_Virus" = 1.0, "Sig_AutoAg" = 1.2, "Nominal" = 0.3), guide="none") +
    scale_alpha_manual(values=c("Sig_Virus" = 0.8, "Sig_AutoAg" = 1.0, "Nominal" = 0.4), guide="none") +
    
    # Nodes
    geom_point(data=nodes_plot, aes(x=x_adj, y=y_adj, size=node_size, fill=node_fill), 
               shape=21, color="grey30", stroke=0.3) +
    scale_fill_identity(name = "Node Type") + 
    scale_size_identity() +
    
    # Labels
    
    # 1. Virus Peptides (Center) - Generally don't label these if too many, but if needed:
    # (Uncomment if you want viral peptide labels, usually they are IDs and crowd the center)
    # geom_text_repel(data=nodes_plot %>% filter(ntype=="Virus"), 
    #                 aes(x=x_adj, y=y_adj, label=name), size=2, segment.size=0.2) +
    
    # 2. Auto-Ag (Middle Ring)
    geom_text_repel(data=nodes_plot %>% filter(ntype=="Auto-Ag"), 
                    aes(x=x_adj, y=y_adj, label=name), 
                    size=4, fontface="bold", bg.color = "white", bg.r = 0.15, 
                    max.overlaps = 50) +
    
    # 3. Disease (Outer Ring)
    geom_text_repel(data=nodes_plot %>% filter(ntype=="Disease"), 
                    aes(x=x_adj, y=y_adj, label=name), 
                    size=3.8, box.padding = 0.5, color='black', max.overlaps = 50) +
    
    theme_void() + 
    theme(legend.position = "bottom", 
          legend.text  = element_text(size = 12),
          legend.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(hjust = 0.5, face="bold", size=16)) +
    labs(title = paste0(target_virus, " network (", outcome_type, ")\nVirus (Center) | Auto-Ag (Middle) | Disease (Outer)"))
  
  safe_virus_name <- gsub("[^A-Za-z0-9]", "_", target_virus)
  filename <- sprintf("network_%s_%s.png", safe_virus_name, tolower(outcome_type))
  ggsave(file.path(opt$out_dir, filename), p, width = 13, height = 13, dpi = 300, bg='white')
  cat(sprintf("  Saved: %s\n", filename))
}



# 4. Core Plotting Function for LAB
# ----------------------
run_network_lab <- function(target_virus, outcome_type, human_lab, virus_lab, sig_ab_full) {

  # --- Step A: Filter Data (NOMINAL P < 0.05) ---
  
  # 1. Virus -> Auto-Ag
  edges_vh <- sig_ab_full %>%
    filter(taxon_species == target_virus) %>%
    filter(!is.na(P), P < 0.05) %>%
    transmute(from = var_id, to = gene_symbol, 
              pval = P, fdr = P.adj, 
              etype = "Virus→Auto-Ag")
  
  if(nrow(edges_vh) == 0) {
    cat(sprintf("  Skipping %s: No nominal Virus-Ab signals.\n", target_virus))
    return(NULL)
  }
  
  relevant_auto_ags <- unique(edges_vh$to)
  
  # 2. Auto-Ag -> LAB
  edges_hd <- human_lab %>%
    filter(gene_symbol %in% relevant_auto_ags) %>%
    filter(!is.na(P_all), P_all < 0.05) %>%
    transmute(from = gene_symbol, to = lab_test, 
              pval = P_all, fdr = p.adj_all, 
              etype = "Auto-Ag→LAB")
  
  # 3. Virus -> LAB
  edges_vd <- virus_lab %>%
    filter(taxon_species == target_virus) %>%
    filter(!is.na(P_all), P_all < 0.05) %>%
    transmute(from = UniProt_acc, to = lab_test, 
              pval = P_all, fdr = p.adj_all, 
              etype = "Virus→LAB")
  
  edges_all <- bind_rows(edges_vh, edges_hd, edges_vd)
  
  if(nrow(edges_all) == 0) return(NULL)
  
  # --- Step B: Edge Classification ---
  
  edges_all <- edges_all %>%
    mutate(
      logFDR = -log10(fdr),
      edge_class = case_when(
        etype == "Auto-Ag→LAB" & fdr < 0.05 ~ "Sig_AutoAg",
        (etype == "Virus→Auto-Ag" | etype == "Virus→LAB") & fdr < 0.05 ~ "Sig_Virus",
        TRUE ~ "Nominal"
      ),
      draw_order = case_when(
        edge_class == "Nominal" ~ 1,
        edge_class == "Sig_Virus" ~ 2,
        edge_class == "Sig_AutoAg" ~ 3
      )
    ) %>%
    arrange(draw_order)
  
  # --- Step C: Node definitions & Sizing ---
  
  node_names <- unique(c(edges_all$from, edges_all$to))
  nodes <- tibble::tibble(name = node_names) %>%
    mutate(
      ntype = case_when(
        name %in% sig_ab_full$var_id ~ "Virus", 
        name %in% sig_ab_full$gene_symbol ~ "Auto-Ag",
        name %in% unique(c(human_lab$lab_test, virus_lab$lab_test)) ~ "LAB",
        TRUE ~ "Other"
      )
    )
  
  # Sizing logic
  ha_pep <- edges_vh %>% count(to, name = "n_pep") %>% rename(name = to)
  ha_dis <- edges_hd %>% count(from, name = "n_dis") %>% rename(name = from)
  ha_size <- full_join(ha_pep, ha_dis, by = "name") %>% 
    mutate(across(c(n_pep, n_dis), ~replace_na(.,0)), ha_total = n_pep + n_dis)
  
  dis_v <- edges_vd %>% count(to, name = "n_from_viral") %>% rename(name = to)
  dis_h <- edges_hd %>% count(to, name = "n_from_human") %>% rename(name = to)
  dis_size <- full_join(dis_v, dis_h, by = "name") %>% 
    mutate(across(c(n_from_viral, n_from_human), ~replace_na(.,0)), dis_total = n_from_viral + n_from_human)
  
  edges_hd.sig = subset(edges_hd, edges_hd$fdr <0.05)
  diseases_with_human <- unique(edges_hd.sig$to)
  edges_vd.sig = subset(edges_vd, edges_vd$fdr <0.05)
  diseases_with_virus <- unique(edges_vd.sig$to)
  
  nodes_final <- nodes %>%
    left_join(ha_size, by="name") %>%
    left_join(dis_size, by="name") %>%
    mutate(
      node_fill = case_when(
        ntype == "Auto-Ag" ~ "skyblue",
        ntype == "Virus" ~ "green",
        ntype == "LAB" & name %in% c(diseases_with_human,diseases_with_virus)  ~ "red",
        ntype == "LAB" ~ "pink",
        TRUE ~ "grey70"
      ),
      node_size = case_when(
        ntype == "Auto-Ag" ~ scales::rescale(replace_na(ha_total, 0), to=c(5, 15)), 
        ntype == "LAB" ~ scales::rescale(replace_na(dis_total, 0), to=c(4, 10)),
        ntype == "Virus"   ~ 2.5,
        TRUE ~ 2.5
      )
    )
  nodes_final$node_size[is.na(nodes_final$node_size)] <- 3
  
  # --- Step D: Concentric Layout Calculation ---
  # CONFIGURATION:
  # 1. Inner Core (R=2): Virus Peptides
  # 2. Middle Ring (R=6): Auto-Ags
  # 3. Outer Ring (R=10): LAB
  
  r_inner  <- 2
  r_middle <- 6
  r_outer  <- 10
  
  nodes_plot <- nodes_final
  
  # D.1: LAB (Outer Ring - Grouped)
  dis_subset <- nodes_plot %>% filter(ntype == "LAB")
  if(nrow(dis_subset) > 0) {
    dis_df <- dis_subset %>% left_join(lab_lookup, by = c("name" = "LAB")) %>%
      mutate(Category = replace_na(Category, "Other")) %>%
      arrange(Category, name)
    
    dis_df$cat_id <- as.numeric(factor(dis_df$Category, levels=unique(dis_df$Category)))
    gap_size <- 2
    dis_df <- dis_df %>%
      mutate(
        row_idx = row_number(),
        total_slots = n() + (max(cat_id) * gap_size),
        pos_idx = row_idx + (cat_id * gap_size),
        theta = (pos_idx / max(pos_idx)) * 2 * pi,
        x_adj = r_outer * cos(theta),
        y_adj = r_outer * sin(theta)
      )
    
    nodes_plot <- nodes_plot %>% 
      left_join(dis_df %>% select(name, x_adj, y_adj), by="name")
  } 
  
  # D.2: Auto-Ags (Middle Ring - Radius 6) -> MOVED FROM CENTER TO MIDDLE
  ag_subset <- nodes_plot %>% filter(ntype == "Auto-Ag")
  if(nrow(ag_subset) > 0) {
    theta_ag <- seq(0, 2*pi, length.out = nrow(ag_subset) + 1)[1:nrow(ag_subset)]
    ag_df <- tibble(
      name = ag_subset$name,
      x_adj = r_middle * cos(theta_ag),
      y_adj = r_middle * sin(theta_ag)
    )
    
    nodes_plot <- nodes_plot %>% 
      rows_update(ag_df, by="name")
  }
  
  # D.3: Virus Peptides (Inner Core - Radius 2) -> MOVED FROM MIDDLE TO CENTER
  vir_subset <- nodes_plot %>% filter(ntype == "Virus")
  if(nrow(vir_subset) > 0) {
    theta_vir <- seq(0, 2*pi, length.out = nrow(vir_subset) + 1)[1:nrow(vir_subset)]
    vir_df <- tibble(
      name = vir_subset$name,
      x_adj = r_inner * cos(theta_vir),
      y_adj = r_inner * sin(theta_vir)
    )
    
    nodes_plot <- nodes_plot %>% 
      rows_update(vir_df, by="name")
  }
  
  # Fallback for unconnected nodes
  nodes_plot$x_adj[is.na(nodes_plot$x_adj)] <- 0
  nodes_plot$y_adj[is.na(nodes_plot$y_adj)] <- 0
  
  # Map coordinates to edges
  edges_plot <- edges_all %>%
    left_join(nodes_plot[, c("name","x_adj","y_adj")], by=c("from"="name")) %>% rename(x_from=x_adj, y_from=y_adj) %>%
    left_join(nodes_plot[, c("name","x_adj","y_adj")], by=c("to"="name")) %>% rename(x_to=x_adj, y_to=y_adj)
  
  # --- Step E: Plotting ---
  set.seed(2026)
  p <- ggplot() +
    # Edges
    geom_segment(data=edges_plot, 
                 aes(x=x_from, y=y_from, xend=x_to, yend=y_to, 
                     color=edge_class, linewidth=edge_class, alpha=edge_class), 
                 lineend = "round") +
    
    # Scales
    scale_color_manual(values=c("Sig_Virus"  = "darkgreen", 
                                "Sig_AutoAg" = "red", 
                                "Nominal"    = "grey50"), name="Significance") +
    
    scale_linewidth_manual(values=c("Sig_Virus" = 1.0, "Sig_AutoAg" = 1.2, "Nominal" = 0.3), guide="none") +
    scale_alpha_manual(values=c("Sig_Virus" = 0.8, "Sig_AutoAg" = 1.0, "Nominal" = 0.4), guide="none") +
    
    # Nodes
    geom_point(data=nodes_plot, aes(x=x_adj, y=y_adj, size=node_size, fill=node_fill), 
               shape=21, color="grey30", stroke=0.3) +
    scale_fill_identity(name = "Node Type") + 
    scale_size_identity() +
    
    # Labels
    # 1. Auto-Ag (Middle Ring)
    geom_text_repel(data=nodes_plot %>% filter(ntype=="Auto-Ag"), 
                    aes(x=x_adj, y=y_adj, label=name), 
                    size=4, fontface="bold", bg.color = "white", bg.r = 0.15, 
                    max.overlaps = 50) +
    
    # 2. LAB (Outer Ring)
    geom_text_repel(data=nodes_plot %>% filter(ntype=="LAB"), 
                    aes(x=x_adj, y=y_adj, label=name), 
                    size=3.8, box.padding = 0.5, color='black', max.overlaps = 50) +
    
    theme_void() + 
    theme(legend.position = "bottom", 
          legend.text  = element_text(size = 12),
          legend.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(hjust = 0.5, face="bold", size=16)) +
    labs(title = paste0(target_virus, " network \nVirus (Center) | Auto-Ag (Middle) | LAB (Outer)"))
  
  safe_virus_name <- gsub("[^A-Za-z0-9]", "_", target_virus)
  filename <- sprintf("network_%s_%s.png", safe_virus_name, tolower(outcome_type))
  ggsave(file.path(opt$out_dir, filename), p, width = 13, height = 13, dpi = 300, bg='white')
  cat(sprintf("  Saved: %s\n", filename))
}

# 5. Execution Loop
# ----------------------
all_viruses <- unique(sig_ab$taxon_species)
cat(sprintf("Found %d virus species to process.\n", length(all_viruses)))

for (virus in all_viruses) {
  try({ run_network_analysis(virus, "Incidence", human_dx.inc, virus_dx.inc, sig_ab) })
  try({ run_network_analysis(virus, "Prevalence", human_dx.pre, virus_dx.pre, sig_ab) })
  # try({ run_network_analysis(virus, "Incidence+Prevalence", human_dx.com, virus_dx.com, sig_ab) })
}

for (virus in all_viruses) {
  try({ run_network_lab(virus, "LAB", human_lab, virus_lab, sig_ab) })
}

cat("Done.\n")
