#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(hrbrthemes)
  library(ggplot2)
  library(ggridges)
  library(viridis)
  source("R/utils_common.R"); source("R/utils_qtl.R")
})

# ---------- CLI ----------
opt_list <- list(
  make_option("--llf_promax", type="character", help="LLF VirSIGHT Promax CSV"),
  make_option("--llf_varscore", type="character", help="LLF VirSIGHT VARscores CSV"),
  make_option("--abc_promax", type="character", help="ABC VirSIGHT Promax CSV"),
  make_option("--abc_varscore", type="character", help="ABC VirSIGHT VARscores CSV"),
  make_option("--leo_promax", type="character", help="LEO VirSIGHT Promax CSV"),
  make_option("--leo_varscore", type="character", help="LEO VirSIGHT VARscores CSV"),
  make_option("--out_dir", type="character", default="results/Figure2",
              help="Output directory [default %default]"),
  make_option("--min_pep_prev", type="double", default=1,
              help="Minimum peptide prevalence (%) to show in ridgeline [default %default]")
)

opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- helper to compute peptide seroprev + plot ----------
make_ridge <- function(promax_df, cohort, outfile_png, min_prev = 1) {
  # Filter herpes species
  herpes_promax   <- promax_df[grep("Human [a-z]+herpesvirus", taxon_species, perl = TRUE)]
  
  # Pivot to long + binary
  # Sample columns:
  sample_cols_pr <- grep("^(AM_|0|1|L)", names(herpes_promax),   value = TRUE)
  
  herpes_promax_long <- herpes_promax |>
    dplyr::select(UniProt_acc, taxon_genus, taxon_species, product, all_of(sample_cols_pr)) |>
    pivot_longer(cols = all_of(sample_cols_pr), names_to = "Sample_ID", values_to = "react_value") |>
    mutate(react_cutoff_1 = ifelse(react_value > 1, 1, 0))
  
  # Peptide-level seroprevalence (what we ridge-plot)
  herpes_promax_sum <- herpes_promax_long |>
    group_by(UniProt_acc, taxon_genus, taxon_species, product) |>
    summarise(virus_peptide_seroprev = round(mean(react_cutoff_1) * 100, 2), .groups = "drop")
  
  # Species naming + order
  herpes_promax_sum$taxon_species <- species_label_map(herpes_promax_sum$taxon_species)
  herpes_promax_sum$taxon_species <- factor(
    herpes_promax_sum$taxon_species,
    levels = c("HHV-8","EBV","HHV-7","HHV-6B","HHV-6A","CMV","VZV","HSV-2","HSV-1")
  )
  
  # Plot
  df_plot <- herpes_promax_sum |> filter(virus_peptide_seroprev > min_prev)
  
  p <- ggplot(df_plot, aes(x = virus_peptide_seroprev, y = taxon_species, fill = after_stat(x))) +
    geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.01) +
    scale_fill_viridis(name = "Prevalence (%)", option = "C", alpha = 0.7) +
    xlim(0, 100) +
    labs(
      title = sprintf("Reactivity to Herpesvirus Peptides in %s", cohort),
      x = "Prevalence (%)", y = NULL
    ) +
    theme_ipsum() +
    theme(
      legend.position = "none",
      panel.spacing   = unit(0.1, "lines"),
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 20),
      axis.text.y     = element_text(size = 14, colour = "black"),
      axis.text.x     = element_text(size = 14, colour = "black"),
      strip.text.x    = element_text(size = 14, colour = "black"),
      axis.title.x    = element_text(size = 16, colour = "black")
    )
  
  ggsave(outfile_png, p, width = 10, height = 5, dpi = 300, bg='white')
}

# ---------- Load inputs ----------
llf_promax <- fread(opt$llf_promax)
abc_promax <- fread(opt$abc_promax)
leo_promax <- fread(opt$leo_promax)

lec_promax = inner_join(abc_promax, leo_promax[,-c(2:5)], by=c('UniProt_acc'))

# ---------- Make ridgelines ----------
make_ridge(
  promax_df   = llf_promax,
  cohort      = "MGBB-LLF",
  outfile_png = file.path(opt$out_dir, "seroprev_herpes_promax_ridgeplot_llf.png"),
  min_prev    = opt$min_pep_prev
)

make_ridge(
  promax_df   = lec_promax,
  cohort      = "MGBB-LEC",
  outfile_png = file.path(opt$out_dir, "seroprev_herpes_promax_ridgeplot_lec.png"),
  min_prev    = opt$min_pep_prev
)

message("Done. Wrote virus ridgeline plots to: ", normalizePath(opt$out_dir))


#################################################################################
# Human antibodies

make_ridge.fl <- function(fl_df, cohort, outfile_png, min_prev = 1) {
  if (!"var_id" %in% names(fl_df)) fl_df[, var_id := paste0("h", .I)]

  # Pivot to long + binary
  # Sample columns:
  sample_cols_pr <- grep("^(AM_|0|1|L)", names(fl_df),   value = TRUE)

  auto_abs_long <- fl_df |>
    dplyr::select(var_id, UniProt_acc, gene_symbol, all_of(sample_cols_pr)) |>
    pivot_longer(cols = all_of(sample_cols_pr), names_to = "Sample_ID", values_to = "react_value") |>
    mutate(react_cutoff_1 = ifelse(react_value > 1, 1, 0))

  # Peptide-level seroprevalence (what we ridge-plot)
  auto_abs_sum <- auto_abs_long |>
    group_by(var_id, gene_symbol) |>
    summarise(fl_seroprev = round(mean(react_cutoff_1) * 100, 2)) |> ungroup() |>
    mutate(
      gene_symbol = fct_rev(as.factor(gene_symbol))
    )
  
  # Add dummy y variable for a single ridge
  auto_abs_sum$dummy <- "Human Antibodies"
  
  # Plot
  df_plot <- auto_abs_sum |> filter(fl_seroprev > min_prev)
  p <- ggplot(df_plot, aes(x = fl_seroprev, y = `dummy`, fill = after_stat(x))) +
    geom_density_ridges_gradient(scale = 3, rel_min_height = 0.00001) +
    scale_fill_viridis(name = "Prevalence (%)", option = "C", alpha = 0.7) +
    xlim(0, 60) +
    labs(
      title = sprintf("Reactivity to Human Autoantigens in %s", cohort),
      x = "Prevalence (%)", y = NULL
    ) +
    theme_ipsum() +
    theme(
      legend.position = "none",
      panel.spacing   = unit(0.1, "lines"),
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 20),
      axis.text.y     = element_text(size = 14, colour = "black"),
      axis.text.x     = element_text(size = 14, colour = "black"),
      strip.text.x    = element_text(size = 14, colour = "black"),
      axis.title.x    = element_text(size = 16, colour = "black")
    )
  
  ggsave(outfile_png, p, width = 10, height = 5, dpi = 300, bg='white')
}

# ---------- Load inputs ----------
llf_fl <- fread(opt$llf_fl)
abc_fl <- fread(opt$abc_fl)
leo_fl <- fread(opt$leo_fl)

ale_fl = inner_join(abc_fl, leo_fl[,-c(1,3:10)], by=c('pep_id'))

# ---------- Make ridgelines ----------
make_ridge.fl(
  fl_df   = llf_fl,
  cohort      = "MGBB-LLF",
  outfile_png = file.path(opt$out_dir, "seroprev_human_fl_ridgeplot_llf.png"),
  min_prev    = opt$min_pep_prev
)

make_ridge.fl(
  fl_df   = ale_fl,
  cohort      = "MGBB-LEC",
  outfile_png = file.path(opt$out_dir, "seroprev_human_fl_ridgeplot_lec.png"),
  min_prev    = opt$min_pep_prev
)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(hrbrthemes)
  library(ggplot2)
  library(ggridges)
  library(viridis)
  source("R/utils_common.R"); source("R/utils_qtl.R")
})

# ---------- CLI ----------
opt_list <- list(
  make_option("--llf_promax", type="character", help="LLF VirSIGHT Promax CSV"),
  make_option("--llf_varscore", type="character", help="LLF VirSIGHT VARscores CSV"),
  make_option("--abc_promax", type="character", help="ABC VirSIGHT Promax CSV"),
  make_option("--abc_varscore", type="character", help="ABC VirSIGHT VARscores CSV"),
  make_option("--leo_promax", type="character", help="LEO VirSIGHT Promax CSV"),
  make_option("--leo_varscore", type="character", help="LEO VirSIGHT VARscores CSV"),
  make_option("--out_dir", type="character", default="results/Figure2",
              help="Output directory [default %default]"),
  make_option("--min_pep_prev", type="double", default=1,
              help="Minimum peptide prevalence (%) to show in ridgeline [default %default]")
)

opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- helper to compute peptide seroprev + plot ----------
make_ridge <- function(promax_df, cohort, outfile_png, min_prev = 1) {
  # Filter herpes species
  herpes_promax   <- promax_df[grep("Human [a-z]+herpesvirus", taxon_species, perl = TRUE)]
  
  # Pivot to long + binary
  # Sample columns:
  sample_cols_pr <- grep("^(AM_|0|1|L)", names(herpes_promax),   value = TRUE)
  
  herpes_promax_long <- herpes_promax |>
    dplyr::select(UniProt_acc, taxon_genus, taxon_species, product, all_of(sample_cols_pr)) |>
    pivot_longer(cols = all_of(sample_cols_pr), names_to = "Sample_ID", values_to = "react_value") |>
    mutate(react_cutoff_1 = ifelse(react_value > 1, 1, 0))
  
  # Peptide-level seroprevalence (what we ridge-plot)
  herpes_promax_sum <- herpes_promax_long |>
    group_by(UniProt_acc, taxon_genus, taxon_species, product) |>
    summarise(virus_peptide_seroprev = round(mean(react_cutoff_1) * 100, 2), .groups = "drop")
  
  # Species naming + order
  herpes_promax_sum$taxon_species <- species_label_map(herpes_promax_sum$taxon_species)
  herpes_promax_sum$taxon_species <- factor(
    herpes_promax_sum$taxon_species,
    levels = c("HHV-8","EBV","HHV-7","HHV-6B","HHV-6A","CMV","VZV","HSV-2","HSV-1")
  )
  
  # Plot
  df_plot <- herpes_promax_sum |> filter(virus_peptide_seroprev > min_prev)
  
  p <- ggplot(df_plot, aes(x = virus_peptide_seroprev, y = taxon_species, fill = after_stat(x))) +
    geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.01) +
    scale_fill_viridis(name = "Prevalence (%)", option = "C", alpha = 0.7) +
    xlim(0, 100) +
    labs(
      title = sprintf("Reactivity to Herpesvirus Peptides in %s", cohort),
      x = "Prevalence (%)", y = NULL
    ) +
    theme_ipsum() +
    theme(
      legend.position = "none",
      panel.spacing   = unit(0.1, "lines"),
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 20),
      axis.text.y     = element_text(size = 14, colour = "black"),
      axis.text.x     = element_text(size = 14, colour = "black"),
      strip.text.x    = element_text(size = 14, colour = "black"),
      axis.title.x    = element_text(size = 16, colour = "black")
    )
  
  ggsave(outfile_png, p, width = 10, height = 5, dpi = 300, bg='white')
}

# ---------- Load inputs ----------
llf_promax <- fread(opt$llf_promax)
abc_promax <- fread(opt$abc_promax)
leo_promax <- fread(opt$leo_promax)

lec_promax = inner_join(abc_promax, leo_promax[,-c(2:5)], by=c('UniProt_acc'))

# ---------- Make ridgelines ----------
make_ridge(
  promax_df   = llf_promax,
  cohort      = "MGBB-LLF",
  outfile_png = file.path(opt$out_dir, "seroprev_herpes_promax_ridgeplot_llf.png"),
  min_prev    = opt$min_pep_prev
)

make_ridge(
  promax_df   = lec_promax,
  cohort      = "MGBB-LEC",
  outfile_png = file.path(opt$out_dir, "seroprev_herpes_promax_ridgeplot_lec.png"),
  min_prev    = opt$min_pep_prev
)

message("Done. Wrote virus ridgeline plots to: ", normalizePath(opt$out_dir))


#################################################################################
# Human antibodies

make_ridge.fl <- function(fl_df, cohort, outfile_png, min_prev = 1) {
  if (!"var_id" %in% names(fl_df)) fl_df[, var_id := paste0("h", .I)]

  # Pivot to long + binary
  # Sample columns:
  sample_cols_pr <- grep("^(AM_|0|1|L)", names(fl_df),   value = TRUE)

  auto_abs_long <- fl_df |>
    dplyr::select(var_id, UniProt_acc, gene_symbol, all_of(sample_cols_pr)) |>
    pivot_longer(cols = all_of(sample_cols_pr), names_to = "Sample_ID", values_to = "react_value") |>
    mutate(react_cutoff_1 = ifelse(react_value > 1, 1, 0))

  # Peptide-level seroprevalence (what we ridge-plot)
  auto_abs_sum <- auto_abs_long |>
    group_by(var_id, gene_symbol) |>
    summarise(fl_seroprev = round(mean(react_cutoff_1) * 100, 2)) |> ungroup() |>
    mutate(
      gene_symbol = fct_rev(as.factor(gene_symbol))
    )
  
  # Add dummy y variable for a single ridge
  auto_abs_sum$dummy <- "Human Antibodies"
  
  # Plot
  df_plot <- auto_abs_sum |> filter(fl_seroprev > min_prev)
  p <- ggplot(df_plot, aes(x = fl_seroprev, y = `dummy`, fill = after_stat(x))) +
    geom_density_ridges_gradient(scale = 3, rel_min_height = 0.00001) +
    scale_fill_viridis(name = "Prevalence (%)", option = "C", alpha = 0.7) +
    xlim(0, 60) +
    labs(
      title = sprintf("Reactivity to Human Autoantigens in %s", cohort),
      x = "Prevalence (%)", y = NULL
    ) +
    theme_ipsum() +
    theme(
      legend.position = "none",
      panel.spacing   = unit(0.1, "lines"),
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 20),
      axis.text.y     = element_text(size = 14, colour = "black"),
      axis.text.x     = element_text(size = 14, colour = "black"),
      strip.text.x    = element_text(size = 14, colour = "black"),
      axis.title.x    = element_text(size = 16, colour = "black")
    )
  
  ggsave(outfile_png, p, width = 10, height = 5, dpi = 300, bg='white')
}

# ---------- Load inputs ----------
llf_fl <- fread(opt$llf_fl)
abc_fl <- fread(opt$abc_fl)
leo_fl <- fread(opt$leo_fl)

ale_fl = inner_join(abc_fl, leo_fl[,-c(1,3:10)], by=c('pep_id'))

# ---------- Make ridgelines ----------
make_ridge.fl(
  fl_df   = llf_fl,
  cohort      = "MGBB-LLF",
  outfile_png = file.path(opt$out_dir, "seroprev_human_fl_ridgeplot_llf.png"),
  min_prev    = opt$min_pep_prev
)

make_ridge.fl(
  fl_df   = ale_fl,
  cohort      = "MGBB-LEC",
  outfile_png = file.path(opt$out_dir, "seroprev_human_fl_ridgeplot_lec.png"),
  min_prev    = opt$min_pep_prev
)

message("Done. Wrote ridgeline plots to: ", normalizePath(opt$out_dir))
