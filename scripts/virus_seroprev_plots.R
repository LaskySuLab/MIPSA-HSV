#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
  library(ggridges)
  library(viridis)
})

# ---------- CLI ----------
opt_list <- list(
  make_option("--llf_promax", type="character", help="LLF VirSIGHT Promax CSV"),
  make_option("--llf_varscore", type="character", help="LLF VirSIGHT VARscores CSV"),
  make_option("--abc_promax", type="character", help="ABC VirSIGHT Promax CSV"),
  make_option("--abc_varscore", type="character", help="ABC VirSIGHT VARscores CSV"),
  make_option("--out_dir", type="character", default="results/Figure2",
              help="Output directory [default %default]"),
  make_option("--min_pep_prev", type="double", default=1,
              help="Minimum peptide prevalence (%) to show in ridgeline [default %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

use_theme_ipsum <- requireNamespace("hrbrthemes", quietly = TRUE)
if (use_theme_ipsum) library(hrbrthemes)

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

# ---------- helper to compute peptide seroprev + plot ----------
make_ridge <- function(promax_df, varscore_df, cohort, outfile_png, min_prev = 1) {
  # Filter herpes species
  herpes_varscore <- varscore_df[grep("Human [a-z]+herpesvirus", taxon_species, perl = TRUE)]
  herpes_promax   <- promax_df[grep("Human [a-z]+herpesvirus", taxon_species, perl = TRUE)]
  
  # Pivot to long + binary
  # Sample columns:
  sample_cols_vs <- grep("^(AM_|10)", names(herpes_varscore), value = TRUE)
  sample_cols_pr <- grep("^(AM_|10)", names(herpes_promax),   value = TRUE)
  
  herpes_varscore_long <- herpes_varscore |>
    dplyr::select(taxon_genus, taxon_species, total_peps, all_of(sample_cols_vs)) |>
    pivot_longer(cols = all_of(sample_cols_vs), names_to = "Sample_ID", values_to = "VARScore") |>
    mutate(
      VARScore_cutoff_0 = ifelse(VARScore > 0, 1, 0),
      VARScore_cutoff_1 = ifelse(VARScore > 1, 1, 0)
    )
  
  herpes_promax_long <- herpes_promax |>
    dplyr::select(UniProt_acc, taxon_genus, taxon_species, product, all_of(sample_cols_pr)) |>
    pivot_longer(cols = all_of(sample_cols_pr), names_to = "Sample_ID", values_to = "react_value") |>
    mutate(react_cutoff_1 = ifelse(react_value > 1, 1, 0))
  
  # Summaries (species-level VARScore—kept for completeness; not plotted)
  hvs_0 <- herpes_varscore_long |>
    group_by(taxon_genus, taxon_species, total_peps) |>
    summarise(virus_seroprev = round(mean(VARScore_cutoff_0)*100, 2), .groups="drop") |>
    mutate(cutoff_label = "Cutoff at 0")
  
  hvs_1 <- herpes_varscore_long |>
    group_by(taxon_genus, taxon_species, total_peps) |>
    summarise(virus_seroprev = round(mean(VARScore_cutoff_1)*100, 2), .groups="drop") |>
    mutate(cutoff_label = "Cutoff at 1")
  
  herpes_varscore_sum <- bind_rows(hvs_0, hvs_1) |>
    mutate(species_label = sprintf("%s (N peptides = %d)", taxon_species, total_peps),
           species_label = forcats::fct_rev(as.factor(species_label)),
           cutoff_label  = forcats::fct_rev(as.factor(cutoff_label)))
  
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
    (if (use_theme_ipsum) hrbrthemes::theme_ipsum() else theme_minimal(base_size = 12)) +
    theme(
      legend.position = "none",
      panel.spacing   = unit(0.1, "lines"),
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 20),
      axis.text.y     = element_text(size = 14, colour = "black"),
      axis.text.x     = element_text(size = 14, colour = "black"),
      strip.text.x    = element_text(size = 14, colour = "black"),
      axis.title.x    = element_text(size = 16, colour = "black")
    )
  
  ggsave(outfile_png, p, width = 10, height = 5, dpi = 300)
}

# ---------- Load inputs ----------
llf_promax <- fread(opt$llf_promax)
llf_vs     <- fread(opt$llf_varscore)
abc_promax <- fread(opt$abc_promax)
abc_vs     <- fread(opt$abc_varscore)

# ---------- Make ridgelines ----------
make_ridge(
  promax_df   = llf_promax,
  varscore_df = llf_vs,
  cohort      = "MGBB-LLF",
  outfile_png = file.path(opt$out_dir, "seroprev_herpes_promax_ridgeplot_llf.png"),
  min_prev    = opt$min_pep_prev
)

make_ridge(
  promax_df   = abc_promax,
  varscore_df = abc_vs,
  cohort      = "MGBB-ABC",
  outfile_png = file.path(opt$out_dir, "seroprev_herpes_promax_ridgeplot_abc.png"),
  min_prev    = opt$min_pep_prev
)

message("Done. Wrote virus ridgeline plots to: ", normalizePath(opt$out_dir))
