#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse); library(data.table); library(dplyr); library(tidyr); library(forcats); library(ggplot2);
  library(ggrepel); library(viridis); library(circlize);  library(ComplexHeatmap); library(grid)
  source("R/utils_common.R"); source("R/utils_qtl.R")
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
  
  # LEC files
  make_option("--lec_all",      type="character", help="LEC GLM all results TSV (lec_hsv_bin_glm_test.tsv)"),
  make_option("--lec_vi_bin",   type="character", help="LEC HSV binary TSV"),
  make_option("--lec_hu_bin",   type="character", help="LEC human FL binary TSV"),
  
  # Replication files
  make_option("--rep_llf",      type="character", help="Replication set from LLF (hsv_bin_fchange_rep_llf.tsv)"),
  make_option("--rep_lec",      type="character", help="Replication rows in LEC (hsv_bin_fchange_rep_lec.tsv)"),
  
  # Figure output dir
  make_option("--out_dir", type="character", default="results/Figure3", help="Output directory [default %default]"),
  
  # Extra plot params
  make_option("--llf_n", type="integer", default=1290, help="LLF N for prevalence denominator [default %default]"),
  make_option("--lec_n", type="integer", default=763,  help="LEC N for prevalence denominator [default %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------- prevalence from binary matrix -------------------------
# - hu_bin: columns var_id + Sample_ID ; vi_bin: columns UniProt_acc + Sample_ID
# - map human var_id -> gene_symbol via hu_fl; virus UniProt_acc -> species via vi_promax

prevalence_human <- function(hu_bin, hu_fl, N) {
  hu_fl$var_id <- if ("var_id" %in% names(hu_fl)) hu_fl$var_id else paste0("h", seq_len(nrow(hu_fl)))
  sample_col <- c("Sample_Id", "Subject_Id")
  mat <- as.data.frame(hu_bin[, setdiff(colnames(hu_bin), sample_col), with = FALSE])
  counts <- colSums(mat, na.rm = TRUE)
  df <- data.frame(fl.gene = names(counts), count = as.numeric(counts))
  df <- df %>% left_join(hu_fl[, c("gene_symbol", "UniProt_acc", "var_id")], by = c("fl.gene" = "var_id"))
  df$gene_symbol <- ifelse(is.na(df$gene_symbol), df$UniProt_acc, df$gene_symbol)
  df$prevalence <- df$count * 100 / N
  df$gene_symbol <- ifelse(
    grepl("^ZNF559-ZNF177,ZNF177", df$gene_symbol), "ZNF177", df$gene_symbol)
  df %>%
    group_by(gene_symbol) %>%
    summarise(Prevalence = max(prevalence, na.rm = TRUE), .groups = "drop")
  }


prevalence_virus <- function(vi_bin, vi_promax, N) {
  sample_col <- c("Sample_Id", "Subject_Id")
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

# diamond data (LLF/LEC variants)
diamond_data_from_rep <- function(rep_dt) {
  rep_dt$gene_symbol <- gsub(" \\s*\\([^\\)]+\\)", "",  rep_dt$gene_symbol)
  rep_dt$gene_symbol <- gsub("GPRASP3,ARMCX5-GPRASP2", "GPRASP3", rep_dt$gene_symbol)
  
  dd <- rep_dt %>%
    mutate(log10_P = -log10(P.adj)) %>%
    group_by(taxon_genus, taxon_species, gene_symbol) %>%
    summarise(
      prevalence = max(Prevalence, na.rm = TRUE),
      Association_Count = n(),
      Max_log10_P = max(log10_P, na.rm = TRUE),
      .groups = "drop"
    )
  dd$taxon_species <- species_label_map(dd$taxon_species)
  dd$taxon_species <- factor(dd$taxon_species, levels = c("HSV-1","HSV-2","VZV","CMV","HHV-6A","HHV-6B","HHV-7","EBV","HHV-8"))
  dd <- dd %>%
    mutate(gene_symbol_label = paste0(gene_symbol, " (",format(round(prevalence, 2), nsmall = 2),")")) %>%
    arrange(desc(Max_log10_P)) %>%
    mutate(gene_symbol_label = factor(gene_symbol_label, levels = unique(gene_symbol_label)))
  dd
}

plot_diamond <- function(dd, title_txt, outfile, cap_size = 100, cap_color = 80, color_label = "-log10(FDR)") {
  p <- ggplot(dd, aes(x = gene_symbol_label, y = taxon_species)) +
    geom_point(aes(size = pmin(Association_Count, cap_size),
                   color = pmin(Max_log10_P, cap_color)),
               shape = 18, alpha = 0.8) +
    scale_size_continuous(range = c(2, 10),
                          breaks = c(1, 3, 5, 10, 30),
                          labels = c("1", "3", "5", "10", "30+")) +
    scale_color_viridis_c(option = "viridis", name = color_label) +
    labs(title = title_txt,
         x = "Gene Symbols of Autoantigens (prevalence %)",
         y = "Herpes viruses",
         size = "Association Count") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
      axis.text.y = element_text(hjust = 1, size = 9),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, size = 13),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      legend.text  = element_text(size = 9),
      legend.title = element_text(size = 9)
    )
  ggsave(outfile, p, width = 13, height = 5, bg='white',dpi = 300)
}

# Manhattan (facet by species, color by Beta, shape by direction)
plot_manhattan <- function(dt_collapse, title_txt, cutoff_p = NULL, lab_thresh, outfile = "manhattan.png") {
  dt_collapse$gene_symbol <- gsub(" \\s*\\([^\\)]+\\)", "",  dt_collapse$gene_symbol)
  
  dt_collapse$taxon_species <- species_label_map(dt_collapse$taxon_species)
  dt_collapse$taxon_species <- factor(dt_collapse$taxon_species, levels = species_levels)
  
  dt_collapse$direction_of_effect <- ifelse(dt_collapse$Beta > 0, "positive", "negative")
  
  p <- ggplot(dt_collapse, aes(x = gene_symbol, y = -log10(P))) +
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
    size = 2, segment.size = 0.1, direction = "both",
    segment.color = 'black', max.overlaps = 15
  )
  
  ggsave(outfile, p, width = 12, height = 5, bg='white', dpi = 300)
}

# Circos / chord plot
plot_circos <- function(rep_summary, prevalence_genes, prevalence_species, title_txt, outfile_png) {
  rep_summary$gene_symbol <- gsub(" \\s*\\([^\\)]+\\)", "",  rep_summary$gene_symbol)
  rep_summary$gene_symbol <- gsub("GPRASP3,ARMCX5-GPRASP2", "GPRASP3", rep_summary$gene_symbol)
  
  chord_data <- rep_summary %>% rename(from = gene_symbol, to = taxon_species)
  chord_data$to <- species_label_map(chord_data$to)
  
  min_val <- round(min(chord_data$log10p, na.rm = TRUE), 2)
  max_val <- round(max(chord_data$log10p, na.rm = TRUE), 2)
  col_fun_chord <- circlize::colorRamp2(
    seq(min_val, max_val, length.out = 7),
    colorRampPalette(c("#FAA491", "#A50F15"))(7)
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
    start.degree = 90,
    gap.after = c(rep(3, length(auto_list)-1), 12, rep(3, length(virus_list)-1), 12),
    circle.margin = 0.3,
    track.margin = c(0.05, 0.05)
  )
  
  set.seed(2026)
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
    title = expression('-log'[10]*'(FDR)'),
    col_fun = col_fun_chord,
    at = seq(min_val, max_val, length.out = 6),
    labels = round(seq(min_val, max_val, length.out = 6), 1),
    legend_gp = gpar(fontsize = 10),
    title_gp  = gpar(fontsize = 11, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    direction = "vertical"
  )
  lgd_discrete <- Legend(
    labels = c("Non-significant"),
    legend_gp = gpar(fill = "grey"),
    title = "Association Type"
  )
  
  draw(packLegend(lgd_continuous, lgd_discrete, direction = "vertical"),
       x = unit(1.0, "npc") - unit(5, "mm"),   # closer to the middle right
       y = unit(0.2, "npc"),                 # vertical center
       just = c("right", "center"))
  
  
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
lec_test <- fread(opt$lec_test)

llf_vi_promax <- fread(opt$llf_vi_promax)
llf_hu_fl     <- fread(opt$llf_hu_fl)
llf_vi_bin <- fread(opt$llf_vi_bin)
llf_hu_bin <- fread(opt$llf_hu_bin)
lec_vi_bin <- fread(opt$lec_vi_bin)
lec_hu_bin <- fread(opt$lec_hu_bin)

rep_llf <- fread(opt$rep_llf)  
rep_lec <- fread(opt$rep_lec)

# ---------- prevalence tables (LLF / LEC) ----------
prev_genes_llf   <- prevalence_human(llf_hu_bin, llf_hu_fl, N = opt$llf_n)
prev_species_llf <- prevalence_virus(llf_vi_bin, llf_vi_promax, N = opt$llf_n)
prev_genes_lec   <- prevalence_human(lec_hu_bin, llf_hu_fl, N = opt$lec_n)
prev_species_lec <- prevalence_virus(lec_vi_bin, llf_vi_promax, N = opt$lec_n)

# ------------------------- DIAMOND (3C / 3E) -------------------------
# FDR P-vale
rep_llf = left_join(rep_llf, prev_genes_llf, by=c('gene_symbol'))

dd_llf <- diamond_data_from_rep(rep_llf)
plot_diamond(
  dd_llf,
  title_txt = "",
  outfile   = file.path(opt$out_dir, "figure3C_diamond_llf_fdr.png"),
  cap_size  = 100, cap_color = 80, color_label = "-log10(FDR)"
)

rep_lec = left_join(rep_lec, prev_genes_llf, by=c('gene_symbol'))
dd_lec <- diamond_data_from_rep(rep_lec)

# reorder LEC to match LLF gene order if overlapping
if (nrow(dd_lec) > 0 && nrow(dd_llf) > 0) {
  dd_lec <- dd_lec %>%
    mutate(key = paste0(taxon_species, gene_symbol)) %>%
    arrange(match(key, paste0(dd_llf$taxon_species, dd_llf$gene_symbol))) %>%
    select(-key)
}

dd_lec$gene_symbol_label = factor(dd_lec$gene_symbol_label, levels = unique(dd_lec$gene_symbol_label))

plot_diamond(
  dd_lec,
  title_txt = "",
  outfile   = file.path(opt$out_dir, "figure3E_diamond_lec_fdr.png"),
  cap_size  = 100, cap_color = 12, color_label = "-log10(FDR)"
)

# ------------------------- CIRCOS (3D / 3F) -------------------------
# LLF circos summary from replicated set
rep_sum_llf <- rep_llf %>%
  group_by(taxon_species, gene_symbol) %>%
  summarise(association_count = n(),
            min_adj_pval = min(P.adj, na.rm = TRUE), .groups = "drop") %>%
  mutate(log10p = -log10(min_adj_pval))
rep_sum_llf$taxon_species <- species_label_map(rep_sum_llf$taxon_species)

prev_genes_llf$gene_symbol <- gsub(" \\s*\\([^\\)]+\\)", "",  prev_genes_llf$gene_symbol)
prev_genes_llf$gene_symbol <- gsub("GPRASP3,ARMCX5-GPRASP2", "GPRASP3",  prev_genes_llf$gene_symbol)
prev_genes_named_llf   <- setNames(prev_genes_llf$Prevalence, prev_genes_llf$gene_symbol)
prev_species_named_llf <- setNames(prev_species_llf$Prevalence, prev_species_llf$taxon_species)

plot_circos(
  rep_summary       = rep_sum_llf,
  prevalence_genes  = prev_genes_named_llf,
  prevalence_species= prev_species_named_llf,
  title_txt         = "",
  outfile_png       = file.path(opt$out_dir, "figure3D_circos_llf_fdr.png")
)

# LEC circos summary from LEC replicated rows
rep_sum_lec <- rep_lec %>%
  group_by(taxon_species, gene_symbol) %>%
  summarise(association_count = n(),
            min_pval = min(P, na.rm = TRUE), .groups = "drop") %>%
  mutate(log10p = -log10(min_pval))
rep_sum_lec$taxon_species <- species_label_map(rep_sum_lec$taxon_species)

prev_genes_lec$gene_symbol <- gsub(" \\s*\\([^\\)]+\\)", "",  prev_genes_lec$gene_symbol)
prev_genes_lec$gene_symbol <- gsub("GPRASP3,ARMCX5-GPRASP2", "GPRASP3", prev_genes_lec$gene_symbol)
prev_genes_named_lec   <- setNames(prev_genes_lec$Prevalence, prev_genes_lec$gene_symbol)
prev_species_named_lec <- setNames(prev_species_lec$Prevalence, prev_species_lec$taxon_species)

plot_circos(
  rep_summary       = rep_sum_lec,
  prevalence_genes  = prev_genes_named_lec,
  prevalence_species= prev_species_named_lec,
  title_txt         = "",
  outfile_png       = file.path(opt$out_dir, "figure3F_circos_lec_fdr.png")
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
llf_collapse1 <- na.omit(llf_collapse)

if ("gene_symbol" %in% names(llf_collapse1)) {
  llf_collapse1$gene_symbol <- ifelse(
    grepl("^TSPY3,TSPY10,LOC", llf_collapse1$gene_symbol),
    "TSPY3, TSPY10",
    llf_collapse1$gene_symbol
  )
}

plot_manhattan(
  llf_collapse1,
  title_txt = "",
  cutoff_p  = llf_cutoff_p, lab_thresh = 1e-10,
  outfile   = file.path(opt$out_dir, "figure3A_manhattan_llf.png")
)

# LEC: restrict to QTLs overlapping LLF, then collapse
lec_res1 <- lec_test

lec_cutoff_p <- {
  sig_idx <- which(lec_res1$P.adj <= 0.05)
  if (length(sig_idx)) max(lec_res1$P[sig_idx], na.rm = TRUE) else NULL
}
lec_collapse <- collapse_minP(lec_res1)
plot_manhattan(
  lec_collapse,
  title_txt = "",
  cutoff_p  = lec_cutoff_p, lab_thresh = lec_cutoff_p,
  outfile   = file.path(opt$out_dir, "figure3B_manhattan_lec.png")
)

message("Figure 3 panels written to: ", normalizePath(opt$out_dir))

