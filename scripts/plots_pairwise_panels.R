#!/usr/bin/env Rscript
# Figure 3 panels: Manhattan (LLF/ABC), Diamond (LLF/ABC), Circos (LLF/ABC)
# Inputs expected from earlier steps (run_hu_vs_hsv_glm.R, build_binary_matrices.R)
#
# Outputs (default to results/Figure3):
#   figure3A_manhattan_llf.png
#   figure3B_manhattan_abc.png
#   figure3C_diamond_llf.png
#   figure3D_circos_llf.png
#   figure3E_diamond_abc.png
#   figure3F_circos_abc.png
#
# Notes:
# - Species names are normalized (HSV-1/2, VZV, CMV, HHV-6A/B, HHV-7, EBV, HHV-8)
# - Diamond plots use replicated set (LLF sig replicated in ABC)
# - Circos uses replicated summary + prevalence labels on sectors

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
  library(ggrepel)
  library(viridis)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
})

# ------------------------- CLI -------------------------
opt_list <- list(
  # LLF files
  make_option("--llf_all",      type="character", help="LLF GLM all results TSV (llf_hsv_bin_glm_all.tsv)"),
  make_option("--llf_sig",      type="character", help="LLF GLM FDR sig TSV (llf_hsv_bin_glm_sig.tsv)"),
  make_option("--llf_vi_promax",type="character", help="LLF VirSIGHT Promax CSV for annotation (for species mapping)"),
  make_option("--llf_hu_fl",    type="character", help="LLF HuSIGHT FL CSV for gene_symbol/var_id (h1..hN)"),
  make_option("--llf_vi_bin",   type="character", help="LLF HSV binary TSV"),
  make_option("--llf_hu_bin",   type="character", help="LLF human FL binary TSV"),
  
  # ABC files
  make_option("--abc_all",      type="character", help="ABC GLM all results TSV (abc_hsv_bin_glm_all.tsv)"),
  make_option("--abc_sig",      type="character", help="ABC GLM FDR sig TSV (abc_hsv_bin_glm_sig.tsv)"),
  make_option("--abc_vi_promax",type="character", help="ABC VirSIGHT Promax CSV for annotation"),
  make_option("--abc_hu_fl",    type="character", help="ABC HuSIGHT FL CSV for gene_symbol/var_id"),
  make_option("--abc_vi_bin",   type="character", help="ABC HSV binary TSV"),
  make_option("--abc_hu_bin",   type="character", help="ABC human FL binary TSV"),
  
  # Replication files
  make_option("--rep_llf",      type="character", help="Replication set from LLF (rep_hsv.sig)"),
  make_option("--rep_abc",      type="character", help="Replication rows in ABC (abc_hsv.rep)"),
  
  # Figure output dir
  make_option("--out_dir", type="character", default="results/Figure3", help="Output directory [default %default]"),
  
  # Extra plot params
  make_option("--llf_n", type="integer", default=1290, help="LLF N for prevalence denominator [default %default]"),
  make_option("--abc_n", type="integer", default=300,  help="ABC N for prevalence denominator [default %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- helpers -------------------------
species_label_map <- function(x) {
  x <- gsub('Human alphaherpesvirus 1', 'HSV-1', x)
  x <- gsub('Human alphaherpesvirus 2', 'HSV-2', x)
  x <- gsub('Human alphaherpesvirus 3', 'VZV', x)
  x <- gsub('Human betaherpesvirus 5', 'CMV', x)
  x <- gsub('Human betaherpesvirus 6A', 'HHV-6A', x)
  x <- gsub('Human betaherpesvirus 6B', 'HHV-6B', x)
  x <- gsub('Human betaherpesvirus 7', 'HHV-7', x)
  x <- gsub('Human gammaherpesvirus 4', 'EBV', x)
  x <- gsub('Human gammaherpesvirus 8', 'HHV-8', x)
  x
}
species_levels <- c("HSV-1","HSV-2","VZV","CMV","HHV-6A","HHV-6B","HHV-7","EBV","HHV-8")

cap_vals <- function(x, cap) pmin(x, cap)

# collapse (antibody x product) to min P for Manhattan
collapse_minP <- function(dt) {
  dt %>%
    mutate(product1 = paste0(antibody, "_", product)) %>%
    group_by(product1) %>%
    slice_min(order_by = P, with_ties = FALSE) %>%
    ungroup()
}

# prevalence from binary matrix
# - hu_bin: columns var_id + Sample_ID ; vi_bin: columns UniProt_acc + Sample_ID
# - map human var_id -> gene_symbol via hu_fl; virus UniProt_acc -> species via vi_promax
prevalence_human <- function(hu_bin, hu_fl, N) {
  hu_fl$var_id <- if ("var_id" %in% names(hu_fl)) hu_fl$var_id else paste0("h", seq_len(nrow(hu_fl)))
  sample_col <- "Sample_ID"
  mat <- as.data.frame(hu_bin[, setdiff(colnames(hu_bin), sample_col), with = FALSE])
  counts <- colSums(mat, na.rm = TRUE)
  df <- data.frame(fl.gene = names(counts), count = as.numeric(counts))
  df <- df %>% left_join(hu_fl[, c("gene_symbol", "UniProt_acc", "var_id")], by = c("fl.gene" = "var_id"))
  df$gene_symbol <- ifelse(is.na(df$gene_symbol), df$UniProt_acc, df$gene_symbol)
  df$prevalence <- df$count * 100 / N
  df %>%
    group_by(gene_symbol) %>%
    summarise(Prevalence = max(prevalence, na.rm = TRUE), .groups = "drop")
}

prevalence_virus <- function(vi_bin, vi_promax, N) {
  sample_col <- "Sample_ID"
  mat <- as.data.frame(vi_bin[, setdiff(colnames(vi_bin), sample_col), with = FALSE])
  counts <- colSums(mat, na.rm = TRUE)
  df <- data.frame(virus_pro = names(counts), count = as.numeric(counts))
  df <- df %>% left_join(vi_promax[, c("UniProt_acc","taxon_genus","taxon_species")], by = c("virus_pro" = "UniProt_acc"))
  df$taxon_species <- species_label_map(df$taxon_species)
  df$prevalence <- df$count * 100 / N
  df %>%
    group_by(taxon_species) %>%
    summarise(Prevalence = max(prevalence, na.rm = TRUE), .groups = "drop")
}

# diamond data (LLF/ABC variants)
diamond_data_from_rep <- function(rep_dt, is_llf = TRUE) {
  dd <- rep_dt %>%
    mutate(log10_P = if (is_llf) -log10(P.adj) else -log10(P)) %>%
    group_by(taxon_genus, taxon_species, gene_symbol, Prevalence) %>%
    summarise(
      Association_Count = n(),
      Max_log10_P = max(log10_P, na.rm = TRUE),
      .groups = "drop"
    )
  dd$taxon_species <- species_label_map(dd$taxon_species)
  dd$taxon_species <- factor(dd$taxon_species, levels = c("HSV-1","HSV-2","VZV","CMV","HHV-6A","HHV-6B","EBV"))
  dd <- dd %>%
    mutate(gene_symbol_label = paste0(gene_symbol, " (", format(round(Prevalence, 2), nsmall = 2), ")")) %>%
    arrange(desc(Max_log10_P)) %>%
    mutate(gene_symbol_label = factor(gene_symbol_label, levels = unique(gene_symbol_label)))
  dd
}

plot_diamond <- function(dd, title_txt, outfile, cap_size = 100, cap_color = 80, color_label = "-log10(adjusted P-val)") {
  p <- ggplot(dd, aes(x = gene_symbol_label, y = taxon_species)) +
    geom_point(aes(size = pmin(Association_Count, cap_size),
                   color = pmin(Max_log10_P, cap_color)),
               shape = 18, alpha = 0.8) +
    scale_size_continuous(range = c(2, 10),
                          breaks = c(1, 3, 5, 10, 30),
                          labels = c("1", "3", "5", "10", "30+")) +
    scale_color_viridis_c(option = "viridis", name = color_label) +
    labs(title = title_txt,
         x = "Gene Symbols of Auto-antibodies (prevalence %)",
         y = "Herpes viruses",
         size = "Association Count") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(hjust = 1, size = 10),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 15),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      legend.text  = element_text(size = 11),
      legend.title = element_text(size = 12)
    )
  ggsave(outfile, p, width = 12, height = 5, dpi = 300)
}

# Manhattan (facet by species, color by Beta, shape by direction)
plot_manhattan <- function(dt_collapse, title_txt, cutoff_p = NULL, lab_thresh, outfile = "manhattan.png") {
  dt_collapse$taxon_species <- species_label_map(dt_collapse$taxon_species)
  dt_collapse$taxon_species <- factor(dt_collapse$taxon_species, levels = species_levels)
  
  dt_collapse$direction_of_effect <- ifelse(dt_collapse$Beta > 0, "positive", "negative")
  
  p <- ggplot(dt_collapse, aes(x = antibody, y = -log10(P))) +
    facet_grid(. ~ taxon_species, scales = "free_x") +
    geom_point(aes(color = Beta, shape = direction_of_effect), size = 2.5, alpha = 1) +
    scale_shape_manual(values = c("negative" = 25, "positive" = 17)) +
    scale_color_gradient2(low = "#2b83ba", mid = "white", high = "#d7191c",
                          midpoint = 0, limits = c(-5, 5), oob = scales::squish) +
    { if (!is.null(cutoff_p)) geom_hline(yintercept = -log10(cutoff_p), color = "black", linetype = "dashed") } +
    geom_hline(yintercept = -log10(0.05), color = "black", linetype = "solid") +
    labs(title = title_txt,
         y = expression('-log'[10]*'('*italic(P)*'-value)'),
         color = "Effect Size",
         shape = "Direction of Effect") +
    theme_bw() +
    theme(
      panel.border  = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.text  = element_text(size = 10),
      plot.title   = element_text(hjust = 0.5, size = 16),
      strip.text.x = element_text(size = 12),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y  = element_text(size = 10),
      axis.title.y = element_text(size = 12)
    )

  p <- p + ggrepel::geom_label_repel(
    data = subset(dt_collapse, P < lab_thresh),
    aes(label = gene_symbol),
    size = 2.0, segment.size = 0.3, direction = "both",
    segment.color = 'black', max.overlaps = 10
  )
  
  ggsave(outfile, p, width = 12, height = 5, dpi = 300)
}

# Circos / chord plot
plot_circos <- function(rep_summary, prevalence_genes, prevalence_species, title_txt, outfile_png) {
  chord_data <- rep_summary %>% rename(from = gene_symbol, to = taxon_species)
  chord_data$to <- species_label_map(chord_data$to)
  
  min_val <- round(min(chord_data$log10p, na.rm = TRUE), 2)
  max_val <- round(max(chord_data$log10p, na.rm = TRUE), 2)
  col_fun_chord <- circlize::colorRamp2(
    seq(min_val, max_val, length.out = 9),
    colorRampPalette(c("#FEE5D9", "#A50F15"))(9)
  )
  
  auto_list  <- unique(chord_data$from)
  virus_list <- unique(chord_data$to)
  all_sectors <- c(auto_list, virus_list)
  
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(outfile_png, width = 2700, height = 2400, res = 300)
  } else {
    # fallback for many Linux builds (also headless)
    png(outfile_png, width = 2700, height = 2400, res = 300, type = "cairo")
  }
  
  par(mar = c(1,1,1,1))
  circos.clear()
  circos.par(
    start.degree = 110,
    gap.after = c(rep(3, length(auto_list)-1), 12, rep(3, length(virus_list)-1), 12),
    circle.margin = 0.3,
    track.margin = c(0.05, 0.05)
  )
  
  set.seed(18)
  chordDiagram(
    x = chord_data[, c("from","to","association_count")],
    order = all_sectors,
    col = col_fun_chord(chord_data$log10p),
    transparency = 0.3,
    directional = 0,
    annotationTrack = "grid",
    link.lwd = pmax(1, chord_data$log10p / max_val * 3)
  )
  
  # Titles
  grid.text(title_txt, x = unit(0.5, "npc"), y = unit(0.98, "npc"),
            gp = gpar(fontsize = 16, fontface = "bold"))
  
  # Sector labels
  circos.track(track.index = 1, panel.fun = function(x, y) {
    sector.index <- get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 1.1, sector.index,
                facing = 'clockwise', niceFacing = TRUE,
                adj = c(0, 0.5), cex = 0.8, font = 2)
  }, bg.border = NA)
  
  # Legends
  lgd_continuous <- Legend(
    title = expression('-log'[10]*'(P)'),
    col_fun = col_fun_chord,
    at = seq(min_val, max_val, length.out = 6),
    labels = round(seq(min_val, max_val, length.out = 6), 1),
    legend_gp = gpar(fontsize = 12),
    title_gp  = gpar(fontsize = 13, fontface = "bold"),
    labels_gp = gpar(fontsize = 11),
    direction = "vertical"
  )
  draw(lgd_continuous, x = unit(0.95, "npc"), y = unit(0.5, "npc"),
       just = c("right","center"))
  
  # Prevalence labels track
  circos.track(track.index = 1, panel.fun = function(x, y) {
    sector.index <- get.cell.meta.data("sector.index")
    prev_val <- ifelse(sector.index %in% names(prevalence_genes),
                       prevalence_genes[sector.index],
                       prevalence_species[sector.index])
    if (!is.na(prev_val)) {
      circos.text(
        CELL_META$xcenter, CELL_META$ylim[1] - 0.1,
        labels = format(round(prev_val, 2), nsmall = 2),
        facing = "reverse.clockwise", niceFacing = TRUE,
        adj = c(0, 0.5), cex = 0.7, font = 3, col = "black"
      )
    }
  }, bg.border = NA)
  
  dev.off()
  circos.clear()
}

# ------------------------- read inputs -------------------------
llf_all <- fread(opt$llf_all)
llf_sig <- fread(opt$llf_sig)
abc_all <- fread(opt$abc_all)
abc_sig <- fread(opt$abc_sig)

llf_vi_promax <- fread(opt$llf_vi_promax)
llf_hu_fl     <- fread(opt$llf_hu_fl)
abc_vi_promax <- fread(opt$abc_vi_promax)
abc_hu_fl     <- fread(opt$abc_hu_fl)

llf_vi_bin <- fread(opt$llf_vi_bin)
llf_hu_bin <- fread(opt$llf_hu_bin)
abc_vi_bin <- fread(opt$abc_vi_bin)
abc_hu_bin <- fread(opt$abc_hu_bin)

rep_llf <- fread(opt$rep_llf)   # rep_hsv.sig (LLF replicated set)
rep_abc <- fread(opt$rep_abc)   # abc_hsv.rep (matching ABC rows)

# ---------- prevalence tables (LLF / ABC) ----------
prev_genes_llf   <- prevalence_human(llf_hu_bin, llf_hu_fl, N = opt$llf_n)
prev_species_llf <- prevalence_virus(llf_vi_bin, llf_vi_promax, N = opt$llf_n)
prev_genes_abc   <- prevalence_human(abc_hu_bin, abc_hu_fl, N = opt$abc_n)
prev_species_abc <- prevalence_virus(abc_vi_bin, abc_vi_promax, N = opt$abc_n)

# ------------------------- DIAMOND (3C / 3E) -------------------------
dd_llf <- diamond_data_from_rep(rep_llf, is_llf = TRUE)
plot_diamond(
  dd_llf,
  title_txt = "Autoantibodies Associated with Herpesviruses in MGBB-LLF",
  outfile   = file.path(opt$out_dir, "figure3C_diamond_llf.png"),
  cap_size  = 100, cap_color = 80, color_label = "-log10(adjusted P-val)"
)

dd_abc <- diamond_data_from_rep(rep_abc, is_llf = FALSE)
# reorder ABC to match LLF gene order if overlapping
if (nrow(dd_abc) > 0 && nrow(dd_llf) > 0) {
  dd_abc <- dd_abc %>%
    mutate(key = paste0(taxon_species, gene_symbol)) %>%
    arrange(match(key, paste0(dd_llf$taxon_species, dd_llf$gene_symbol))) %>%
    select(-key)
}

dd_abc$gene_symbol_label = factor(dd_abc$gene_symbol_label, levels = unique(dd_abc$gene_symbol_label))

plot_diamond(
  dd_abc,
  title_txt = "Significant Autoantibodies Associated with Herpesviruses in MGBB-ABC",
  outfile   = file.path(opt$out_dir, "figure3E_diamond_abc.png"),
  cap_size  = 100, cap_color = 12, color_label = "-log10(P-val)"
)

# ------------------------- CIRCOS (3D / 3F) -------------------------
# LLF circos summary from replicated set
rep_sum_llf <- rep_llf %>%
  group_by(taxon_species, gene_symbol) %>%
  summarise(association_count = n(),
            min_adj_pval = min(P.adj, na.rm = TRUE), .groups = "drop") %>%
  mutate(log10p = -log10(min_adj_pval))
rep_sum_llf$taxon_species <- species_label_map(rep_sum_llf$taxon_species)

prev_genes_named_llf   <- setNames(prev_genes_llf$Prevalence, prev_genes_llf$gene_symbol)
prev_species_named_llf <- setNames(prev_species_llf$Prevalence, prev_species_llf$taxon_species)

plot_circos(
  rep_summary       = rep_sum_llf,
  prevalence_genes  = prev_genes_named_llf,
  prevalence_species= prev_species_named_llf,
  title_txt         = "LLF: Herpesvirus–Autoantibody QTLs (Replicated)",
  outfile_png       = file.path(opt$out_dir, "figure3D_circos_llf.png")
)

# ABC circos summary from ABC replicated rows
rep_sum_abc <- rep_abc %>%
  group_by(taxon_species, gene_symbol) %>%
  summarise(association_count = n(),
            min_pval = min(P, na.rm = TRUE), .groups = "drop") %>%
  mutate(log10p = -log10(min_pval))
rep_sum_abc$taxon_species <- species_label_map(rep_sum_abc$taxon_species)

prev_genes_named_abc   <- setNames(prev_genes_abc$Prevalence, prev_genes_abc$gene_symbol)
prev_species_named_abc <- setNames(prev_species_abc$Prevalence, prev_species_abc$taxon_species)

plot_circos(
  rep_summary       = rep_sum_abc,
  prevalence_genes  = prev_genes_named_abc,
  prevalence_species= prev_species_named_abc,
  title_txt         = "ABC: Herpesvirus–Autoantibody QTLs (Replicated)",
  outfile_png       = file.path(opt$out_dir, "figure3F_circos_abc.png")
)

# ------------------------- MANHATTAN (3A / 3B) -------------------------
# LLF: collapse by antibody x product, draw lines at FDR-cutoff p and 0.05
llf_cutoff_p <- {
  if (nrow(llf_sig) > 0) {
    sig_idx <- which(llf_all$P.adj <= 0.05)
    if (length(sig_idx)) max(llf_all$P[sig_idx], na.rm = TRUE) else NULL
  } else NULL
}
llf_collapse <- collapse_minP(llf_all)
plot_manhattan(
  llf_collapse,
  title_txt = "Auto-antibodies and Herpesviruses in MGBB-LLF cohort",
  cutoff_p  = llf_cutoff_p, lab_thresh = 1e-15,
  outfile   = file.path(opt$out_dir, "figure3A_manhattan_llf.png")
)

# ABC: restrict to QTLs overlapping LLF, then collapse
abc_all$qtl <- paste0(abc_all$antibody, "_", abc_all$var_id)
llf_all$qtl <- paste0(llf_all$antibody, "_", llf_all$var_id)
overlap_keys <- intersect(abc_all$qtl, llf_all$qtl[llf_all$P.adj <= 0.05])

abc_res1 <- abc_all %>% filter(qtl %in% overlap_keys)
abc_res1$P.adj1 <- p.adjust(abc_res1$P, method = 'fdr')
abc_cutoff_p <- {
  sig_idx <- which(abc_res1$P.adj1 <= 0.05)
  if (length(sig_idx)) max(abc_res1$P[sig_idx], na.rm = TRUE) else NULL
}
abc_collapse <- collapse_minP(abc_res1)
plot_manhattan(
  abc_collapse,
  title_txt = "Auto-antibodies and Herpesviruses in MGBB-ABC cohort",
  cutoff_p  = abc_cutoff_p, lab_thresh = abc_cutoff_p,
  outfile   = file.path(opt$out_dir, "figure3B_manhattan_abc.png")
)

message("Figure 3 panels written to: ", normalizePath(opt$out_dir))

#!/usr/bin/env Rscript
# Figure 3 panels: Manhattan (LLF/ABC), Diamond (LLF/ABC), Circos (LLF/ABC)
# Inputs expected from earlier steps (run_hu_vs_hsv_glm.R, build_binary_matrices.R)
#
# Outputs (default to results/Figure3):
#   figure3A_manhattan_llf.png
#   figure3B_manhattan_abc.png
#   figure3C_diamond_llf.png
#   figure3D_circos_llf.png
#   figure3E_diamond_abc.png
#   figure3F_circos_abc.png
#
# Notes:
# - Species names are normalized (HSV-1/2, VZV, CMV, HHV-6A/B, HHV-7, EBV, HHV-8)
# - Diamond plots use replicated set (LLF sig replicated in ABC)
# - Circos uses replicated summary + prevalence labels on sectors

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
  library(ggrepel)
  library(viridis)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
})

# ------------------------- CLI -------------------------
opt_list <- list(
  # LLF files
  make_option("--llf_all",      type="character", help="LLF GLM all results TSV (llf_hsv_bin_glm_all.tsv)"),
  make_option("--llf_sig",      type="character", help="LLF GLM FDR sig TSV (llf_hsv_bin_glm_sig.tsv)"),
  make_option("--llf_vi_promax",type="character", help="LLF VirSIGHT Promax CSV for annotation (for species mapping)"),
  make_option("--llf_hu_fl",    type="character", help="LLF HuSIGHT FL CSV for gene_symbol/var_id (h1..hN)"),
  make_option("--llf_vi_bin",   type="character", help="LLF HSV binary TSV"),
  make_option("--llf_hu_bin",   type="character", help="LLF human FL binary TSV"),
  
  # ABC files
  make_option("--abc_all",      type="character", help="ABC GLM all results TSV (abc_hsv_bin_glm_all.tsv)"),
  make_option("--abc_sig",      type="character", help="ABC GLM FDR sig TSV (abc_hsv_bin_glm_sig.tsv)"),
  make_option("--abc_vi_promax",type="character", help="ABC VirSIGHT Promax CSV for annotation"),
  make_option("--abc_hu_fl",    type="character", help="ABC HuSIGHT FL CSV for gene_symbol/var_id"),
  make_option("--abc_vi_bin",   type="character", help="ABC HSV binary TSV"),
  make_option("--abc_hu_bin",   type="character", help="ABC human FL binary TSV"),
  
  # Replication files
  make_option("--rep_llf",      type="character", help="Replication set from LLF (rep_hsv.sig)"),
  make_option("--rep_abc",      type="character", help="Replication rows in ABC (abc_hsv.rep)"),
  
  # Figure output dir
  make_option("--out_dir", type="character", default="results/Figure3", help="Output directory [default %default]"),
  
  # Extra plot params
  make_option("--llf_n", type="integer", default=1290, help="LLF N for prevalence denominator [default %default]"),
  make_option("--abc_n", type="integer", default=300,  help="ABC N for prevalence denominator [default %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- helpers -------------------------
species_label_map <- function(x) {
  x <- gsub('Human alphaherpesvirus 1', 'HSV-1', x)
  x <- gsub('Human alphaherpesvirus 2', 'HSV-2', x)
  x <- gsub('Human alphaherpesvirus 3', 'VZV', x)
  x <- gsub('Human betaherpesvirus 5', 'CMV', x)
  x <- gsub('Human betaherpesvirus 6A', 'HHV-6A', x)
  x <- gsub('Human betaherpesvirus 6B', 'HHV-6B', x)
  x <- gsub('Human betaherpesvirus 7', 'HHV-7', x)
  x <- gsub('Human gammaherpesvirus 4', 'EBV', x)
  x <- gsub('Human gammaherpesvirus 8', 'HHV-8', x)
  x
}
species_levels <- c("HSV-1","HSV-2","VZV","CMV","HHV-6A","HHV-6B","HHV-7","EBV","HHV-8")

cap_vals <- function(x, cap) pmin(x, cap)

# collapse (antibody x product) to min P for Manhattan
collapse_minP <- function(dt) {
  dt %>%
    mutate(product1 = paste0(antibody, "_", product)) %>%
    group_by(product1) %>%
    slice_min(order_by = P, with_ties = FALSE) %>%
    ungroup()
}

# prevalence from binary matrix
# - hu_bin: columns var_id + Sample_ID ; vi_bin: columns UniProt_acc + Sample_ID
# - map human var_id -> gene_symbol via hu_fl; virus UniProt_acc -> species via vi_promax
prevalence_human <- function(hu_bin, hu_fl, N) {
  hu_fl$var_id <- if ("var_id" %in% names(hu_fl)) hu_fl$var_id else paste0("h", seq_len(nrow(hu_fl)))
  sample_col <- "Sample_ID"
  mat <- as.data.frame(hu_bin[, setdiff(colnames(hu_bin), sample_col), with = FALSE])
  counts <- colSums(mat, na.rm = TRUE)
  df <- data.frame(fl.gene = names(counts), count = as.numeric(counts))
  df <- df %>% left_join(hu_fl[, c("gene_symbol", "UniProt_acc", "var_id")], by = c("fl.gene" = "var_id"))
  df$gene_symbol <- ifelse(is.na(df$gene_symbol), df$UniProt_acc, df$gene_symbol)
  df$prevalence <- df$count * 100 / N
  df %>%
    group_by(gene_symbol) %>%
    summarise(Prevalence = max(prevalence, na.rm = TRUE), .groups = "drop")
}

prevalence_virus <- function(vi_bin, vi_promax, N) {
  sample_col <- "Sample_ID"
  mat <- as.data.frame(vi_bin[, setdiff(colnames(vi_bin), sample_col), with = FALSE])
  counts <- colSums(mat, na.rm = TRUE)
  df <- data.frame(virus_pro = names(counts), count = as.numeric(counts))
  df <- df %>% left_join(vi_promax[, c("UniProt_acc","taxon_genus","taxon_species")], by = c("virus_pro" = "UniProt_acc"))
  df$taxon_species <- species_label_map(df$taxon_species)
  df$prevalence <- df$count * 100 / N
  df %>%
    group_by(taxon_species) %>%
    summarise(Prevalence = max(prevalence, na.rm = TRUE), .groups = "drop")
}

# diamond data (LLF/ABC variants)
diamond_data_from_rep <- function(rep_dt, is_llf = TRUE) {
  dd <- rep_dt %>%
    mutate(log10_P = if (is_llf) -log10(P.adj) else -log10(P)) %>%
    group_by(taxon_genus, taxon_species, gene_symbol, Prevalence) %>%
    summarise(
      Association_Count = n(),
      Max_log10_P = max(log10_P, na.rm = TRUE),
      .groups = "drop"
    )
  dd$taxon_species <- species_label_map(dd$taxon_species)
  dd$taxon_species <- factor(dd$taxon_species, levels = c("HSV-1","HSV-2","VZV","CMV","HHV-6A","HHV-6B","EBV"))
  dd <- dd %>%
    mutate(gene_symbol_label = paste0(gene_symbol, " (", format(round(Prevalence, 2), nsmall = 2), ")")) %>%
    arrange(desc(Max_log10_P)) %>%
    mutate(gene_symbol_label = factor(gene_symbol_label, levels = unique(gene_symbol_label)))
  dd
}

plot_diamond <- function(dd, title_txt, outfile, cap_size = 100, cap_color = 80, color_label = "-log10(adjusted P-val)") {
  p <- ggplot(dd, aes(x = gene_symbol_label, y = taxon_species)) +
    geom_point(aes(size = pmin(Association_Count, cap_size),
                   color = pmin(Max_log10_P, cap_color)),
               shape = 18, alpha = 0.8) +
    scale_size_continuous(range = c(2, 10),
                          breaks = c(1, 3, 5, 10, 30),
                          labels = c("1", "3", "5", "10", "30+")) +
    scale_color_viridis_c(option = "viridis", name = color_label) +
    labs(title = title_txt,
         x = "Gene Symbols of Auto-antibodies (prevalence %)",
         y = "Herpes viruses",
         size = "Association Count") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(hjust = 1, size = 10),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 15),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      legend.text  = element_text(size = 11),
      legend.title = element_text(size = 12)
    )
  ggsave(outfile, p, width = 12, height = 5, dpi = 300)
}

# Manhattan (facet by species, color by Beta, shape by direction)
plot_manhattan <- function(dt_collapse, title_txt, cutoff_p = NULL, lab_thresh, outfile = "manhattan.png") {
  dt_collapse$taxon_species <- species_label_map(dt_collapse$taxon_species)
  dt_collapse$taxon_species <- factor(dt_collapse$taxon_species, levels = species_levels)
  
  dt_collapse$direction_of_effect <- ifelse(dt_collapse$Beta > 0, "positive", "negative")
  
  p <- ggplot(dt_collapse, aes(x = antibody, y = -log10(P))) +
    facet_grid(. ~ taxon_species, scales = "free_x") +
    geom_point(aes(color = Beta, shape = direction_of_effect), size = 2.5, alpha = 1) +
    scale_shape_manual(values = c("negative" = 25, "positive" = 17)) +
    scale_color_gradient2(low = "#2b83ba", mid = "white", high = "#d7191c",
                          midpoint = 0, limits = c(-5, 5), oob = scales::squish) +
    { if (!is.null(cutoff_p)) geom_hline(yintercept = -log10(cutoff_p), color = "black", linetype = "dashed") } +
    geom_hline(yintercept = -log10(0.05), color = "black", linetype = "solid") +
    labs(title = title_txt,
         y = expression('-log'[10]*'('*italic(P)*'-value)'),
         color = "Effect Size",
         shape = "Direction of Effect") +
    theme_bw() +
    theme(
      panel.border  = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.text  = element_text(size = 10),
      plot.title   = element_text(hjust = 0.5, size = 16),
      strip.text.x = element_text(size = 12),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y  = element_text(size = 10),
      axis.title.y = element_text(size = 12)
    )

  p <- p + ggrepel::geom_label_repel(
    data = subset(dt_collapse, P < lab_thresh),
    aes(label = gene_symbol),
    size = 2.0, segment.size = 0.3, direction = "both",
    segment.color = 'black', max.overlaps = 10
  )
  
  ggsave(outfile, p, width = 12, height = 5, dpi = 300)
}

# Circos / chord plot
plot_circos <- function(rep_summary, prevalence_genes, prevalence_species, title_txt, outfile_png) {
  chord_data <- rep_summary %>% rename(from = gene_symbol, to = taxon_species)
  chord_data$to <- species_label_map(chord_data$to)
  
  min_val <- round(min(chord_data$log10p, na.rm = TRUE), 2)
  max_val <- round(max(chord_data$log10p, na.rm = TRUE), 2)
  col_fun_chord <- circlize::colorRamp2(
    seq(min_val, max_val, length.out = 9),
    colorRampPalette(c("#FEE5D9", "#A50F15"))(9)
  )
  
  auto_list  <- unique(chord_data$from)
  virus_list <- unique(chord_data$to)
  all_sectors <- c(auto_list, virus_list)
  
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(outfile_png, width = 2700, height = 2400, res = 300)
  } else {
    # fallback for many Linux builds (also headless)
    png(outfile_png, width = 2700, height = 2400, res = 300, type = "cairo")
  }
  
  par(mar = c(1,1,1,1))
  circos.clear()
  circos.par(
    start.degree = 110,
    gap.after = c(rep(3, length(auto_list)-1), 12, rep(3, length(virus_list)-1), 12),
    circle.margin = 0.3,
    track.margin = c(0.05, 0.05)
  )
  
  set.seed(18)
  chordDiagram(
    x = chord_data[, c("from","to","association_count")],
    order = all_sectors,
    col = col_fun_chord(chord_data$log10p),
    transparency = 0.3,
    directional = 0,
    annotationTrack = "grid",
    link.lwd = pmax(1, chord_data$log10p / max_val * 3)
  )
  
  # Titles
  grid.text(title_txt, x = unit(0.5, "npc"), y = unit(0.98, "npc"),
            gp = gpar(fontsize = 16, fontface = "bold"))
  
  # Sector labels
  circos.track(track.index = 1, panel.fun = function(x, y) {
    sector.index <- get.cell.meta.data("sector.index")
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 1.1, sector.index,
                facing = 'clockwise', niceFacing = TRUE,
                adj = c(0, 0.5), cex = 0.8, font = 2)
  }, bg.border = NA)
  
  # Legends
  lgd_continuous <- Legend(
    title = expression('-log'[10]*'(P)'),
    col_fun = col_fun_chord,
    at = seq(min_val, max_val, length.out = 6),
    labels = round(seq(min_val, max_val, length.out = 6), 1),
    legend_gp = gpar(fontsize = 12),
    title_gp  = gpar(fontsize = 13, fontface = "bold"),
    labels_gp = gpar(fontsize = 11),
    direction = "vertical"
  )
  draw(lgd_continuous, x = unit(0.95, "npc"), y = unit(0.5, "npc"),
       just = c("right","center"))
  
  # Prevalence labels track
  circos.track(track.index = 1, panel.fun = function(x, y) {
    sector.index <- get.cell.meta.data("sector.index")
    prev_val <- ifelse(sector.index %in% names(prevalence_genes),
                       prevalence_genes[sector.index],
                       prevalence_species[sector.index])
    if (!is.na(prev_val)) {
      circos.text(
        CELL_META$xcenter, CELL_META$ylim[1] - 0.1,
        labels = format(round(prev_val, 2), nsmall = 2),
        facing = "reverse.clockwise", niceFacing = TRUE,
        adj = c(0, 0.5), cex = 0.7, font = 3, col = "black"
      )
    }
  }, bg.border = NA)
  
  dev.off()
  circos.clear()
}

# ------------------------- read inputs -------------------------
llf_all <- fread(opt$llf_all)
llf_sig <- fread(opt$llf_sig)
abc_all <- fread(opt$abc_all)
abc_sig <- fread(opt$abc_sig)

llf_vi_promax <- fread(opt$llf_vi_promax)
llf_hu_fl     <- fread(opt$llf_hu_fl)
abc_vi_promax <- fread(opt$abc_vi_promax)
abc_hu_fl     <- fread(opt$abc_hu_fl)

llf_vi_bin <- fread(opt$llf_vi_bin)
llf_hu_bin <- fread(opt$llf_hu_bin)
abc_vi_bin <- fread(opt$abc_vi_bin)
abc_hu_bin <- fread(opt$abc_hu_bin)

rep_llf <- fread(opt$rep_llf)   # rep_hsv.sig (LLF replicated set)
rep_abc <- fread(opt$rep_abc)   # abc_hsv.rep (matching ABC rows)

# ---------- prevalence tables (LLF / ABC) ----------
prev_genes_llf   <- prevalence_human(llf_hu_bin, llf_hu_fl, N = opt$llf_n)
prev_species_llf <- prevalence_virus(llf_vi_bin, llf_vi_promax, N = opt$llf_n)
prev_genes_abc   <- prevalence_human(abc_hu_bin, abc_hu_fl, N = opt$abc_n)
prev_species_abc <- prevalence_virus(abc_vi_bin, abc_vi_promax, N = opt$abc_n)

# ------------------------- DIAMOND (3C / 3E) -------------------------
dd_llf <- diamond_data_from_rep(rep_llf, is_llf = TRUE)
plot_diamond(
  dd_llf,
  title_txt = "Autoantibodies Associated with Herpesviruses in MGBB-LLF",
  outfile   = file.path(opt$out_dir, "figure3C_diamond_llf.png"),
  cap_size  = 100, cap_color = 80, color_label = "-log10(adjusted P-val)"
)

dd_abc <- diamond_data_from_rep(rep_abc, is_llf = FALSE)
# reorder ABC to match LLF gene order if overlapping
if (nrow(dd_abc) > 0 && nrow(dd_llf) > 0) {
  dd_abc <- dd_abc %>%
    mutate(key = paste0(taxon_species, gene_symbol)) %>%
    arrange(match(key, paste0(dd_llf$taxon_species, dd_llf$gene_symbol))) %>%
    select(-key)
}

dd_abc$gene_symbol_label = factor(dd_abc$gene_symbol_label, levels = unique(dd_abc$gene_symbol_label))

plot_diamond(
  dd_abc,
  title_txt = "Significant Autoantibodies Associated with Herpesviruses in MGBB-ABC",
  outfile   = file.path(opt$out_dir, "figure3E_diamond_abc.png"),
  cap_size  = 100, cap_color = 12, color_label = "-log10(P-val)"
)

# ------------------------- CIRCOS (3D / 3F) -------------------------
# LLF circos summary from replicated set
rep_sum_llf <- rep_llf %>%
  group_by(taxon_species, gene_symbol) %>%
  summarise(association_count = n(),
            min_adj_pval = min(P.adj, na.rm = TRUE), .groups = "drop") %>%
  mutate(log10p = -log10(min_adj_pval))
rep_sum_llf$taxon_species <- species_label_map(rep_sum_llf$taxon_species)

prev_genes_named_llf   <- setNames(prev_genes_llf$Prevalence, prev_genes_llf$gene_symbol)
prev_species_named_llf <- setNames(prev_species_llf$Prevalence, prev_species_llf$taxon_species)

plot_circos(
  rep_summary       = rep_sum_llf,
  prevalence_genes  = prev_genes_named_llf,
  prevalence_species= prev_species_named_llf,
  title_txt         = "LLF: Herpesvirus–Autoantibody QTLs (Replicated)",
  outfile_png       = file.path(opt$out_dir, "figure3D_circos_llf.png")
)

# ABC circos summary from ABC replicated rows
rep_sum_abc <- rep_abc %>%
  group_by(taxon_species, gene_symbol) %>%
  summarise(association_count = n(),
            min_pval = min(P, na.rm = TRUE), .groups = "drop") %>%
  mutate(log10p = -log10(min_pval))
rep_sum_abc$taxon_species <- species_label_map(rep_sum_abc$taxon_species)

prev_genes_named_abc   <- setNames(prev_genes_abc$Prevalence, prev_genes_abc$gene_symbol)
prev_species_named_abc <- setNames(prev_species_abc$Prevalence, prev_species_abc$taxon_species)

plot_circos(
  rep_summary       = rep_sum_abc,
  prevalence_genes  = prev_genes_named_abc,
  prevalence_species= prev_species_named_abc,
  title_txt         = "ABC: Herpesvirus–Autoantibody QTLs (Replicated)",
  outfile_png       = file.path(opt$out_dir, "figure3F_circos_abc.png")
)

# ------------------------- MANHATTAN (3A / 3B) -------------------------
# LLF: collapse by antibody x product, draw lines at FDR-cutoff p and 0.05
llf_cutoff_p <- {
  if (nrow(llf_sig) > 0) {
    sig_idx <- which(llf_all$P.adj <= 0.05)
    if (length(sig_idx)) max(llf_all$P[sig_idx], na.rm = TRUE) else NULL
  } else NULL
}
llf_collapse <- collapse_minP(llf_all)
plot_manhattan(
  llf_collapse,
  title_txt = "Auto-antibodies and Herpesviruses in MGBB-LLF cohort",
  cutoff_p  = llf_cutoff_p, lab_thresh = 1e-15,
  outfile   = file.path(opt$out_dir, "figure3A_manhattan_llf.png")
)

# ABC: restrict to QTLs overlapping LLF, then collapse
abc_all$qtl <- paste0(abc_all$antibody, "_", abc_all$var_id)
llf_all$qtl <- paste0(llf_all$antibody, "_", llf_all$var_id)
overlap_keys <- intersect(abc_all$qtl, llf_all$qtl[llf_all$P.adj <= 0.05])

abc_res1 <- abc_all %>% filter(qtl %in% overlap_keys)
abc_res1$P.adj1 <- p.adjust(abc_res1$P, method = 'fdr')
abc_cutoff_p <- {
  sig_idx <- which(abc_res1$P.adj1 <= 0.05)
  if (length(sig_idx)) max(abc_res1$P[sig_idx], na.rm = TRUE) else NULL
}
abc_collapse <- collapse_minP(abc_res1)
plot_manhattan(
  abc_collapse,
  title_txt = "Auto-antibodies and Herpesviruses in MGBB-ABC cohort",
  cutoff_p  = abc_cutoff_p, lab_thresh = abc_cutoff_p,
  outfile   = file.path(opt$out_dir, "figure3B_manhattan_abc.png")
)

message("Figure 3 panels written to: ", normalizePath(opt$out_dir))

