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
  # hrbrthemes optional
})

# ---------- CLI ----------
opt_list <- list(
  make_option("--llf_hu_fl", type="character",
              help="LLF HuSIGHT FullLength (Hits_Fold-Over-Background) CSV"),
  make_option("--abc_hu_fl", type="character",
              help="ABC HuSIGHT FullLength (Hits_Fold-Over-Background) CSV"),
  make_option("--rep_sig", type="character", default=NULL,
              help="Optional: rep_hsv.sig to subset antibodies of interest (not used in plot if NULL)"),
  make_option("--out_dir", type="character", default="results/Figure2",
              help="Output directory [default %default]"),
  make_option("--min_prev", type="double", default=1,
              help="Minimum FL prevalence (%) to show in ridgeline [default %default]"),
  make_option("--xmax", type="double", default=60,
              help="X-axis maximum for prevalence (%) [default %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

use_theme_ipsum <- requireNamespace("hrbrthemes", quietly = TRUE)
if (use_theme_ipsum) library(hrbrthemes)

# ---------- helpers ----------
compute_hu_seroprev <- function(hu_df, sample_prefix) {
  hu_df$var_id <- paste0("h", seq_len(nrow(hu_df)))
  sample_cols <- grep(paste0("^", sample_prefix), names(hu_df), value = TRUE)

  hu_long <- hu_df |>
    dplyr::select(var_id, UniProt_acc, gene_symbol, all_of(sample_cols)) |>
    pivot_longer(cols = all_of(sample_cols), names_to = "Sample_ID", values_to = "react_value") |>
    mutate(react_cutoff_1 = ifelse(react_value > 1, 1, 0))

  hu_sum <- hu_long |>
    group_by(var_id, gene_symbol) |>
    summarise(fl_seroprev = round(mean(react_cutoff_1) * 100, 2), .groups = "drop") |>
    mutate(gene_symbol = fct_rev(as.factor(gene_symbol)))

  hu_sum
}

plot_single_strip <- function(df_sum, title_txt, outfile_png, min_prev = 1, xmax = 60) {
  df_sum$dummy <- "Human Antibodies"

  p <- df_sum |>
    filter(fl_seroprev > min_prev) |>
    ggplot(aes(x = fl_seroprev, y = dummy, fill = after_stat(x))) +
    geom_density_ridges_gradient(scale = 3, rel_min_height = 1e-5) +
    scale_fill_viridis(name = "Prevalence (%)", option = "C", alpha = 0.7) +
    xlim(0, xmax) +
    labs(title = title_txt, x = "Prevalence (%)", y = NULL) +
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
llf_hu <- fread(opt$llf_hu_fl)
abc_hu <- fread(opt$abc_hu_fl)

# Optional subset set (not used in plots unless you want me to change it)
if (!is.null(opt$rep_sig) && file.exists(opt$rep_sig)) {
  rep_sig <- fread(opt$rep_sig)
  rep_abs <- unique(rep_sig$antibody)
} else {
  rep_abs <- NULL
}

# ---------- Compute seroprevalence ----------
# LLF: columns start with "AM_"
llf_sum <- compute_hu_seroprev(llf_hu, sample_prefix = "AM_")
# ABC: columns start with "10"
abc_sum <- compute_hu_seroprev(abc_hu, sample_prefix = "10")

# If you want to constrain to significant antibodies only, uncomment:
# if (!is.null(rep_abs)) {
#   llf_sum <- dplyr::filter(llf_sum, var_id %in% rep_abs)
#   abc_sum <- dplyr::filter(abc_sum, var_id %in% rep_abs)
# }

# ---------- Plot ----------
plot_single_strip(
  df_sum     = llf_sum,
  title_txt  = "Reactivity to Human Autoantibodies in MGBB-LLF",
  outfile_png= file.path(opt$out_dir, "seroprev_human_fl_ridgeplot_asthma.png"),
  min_prev   = opt$min_prev,
  xmax       = opt$xmax
)

plot_single_strip(
  df_sum     = abc_sum,
  title_txt  = "Reactivity to Human Autoantibodies in MGBB-ABC",
  outfile_png= file.path(opt$out_dir, "seroprev_human_fl_ridgeplot_abc.png"),
  min_prev   = opt$min_prev,
  xmax       = opt$xmax
)

message("Done. Wrote human ridgeline plots to: ", normalizePath(opt$out_dir))
