#!/usr/bin/env Rscript
# Build binary matrices for HSV peptides and Human Full-Length antibodies
# Description:
#   - Converts raw FoB data (>1 → 1, else 0) to binary matrices
#   - Removes control samples
#   - Keeps features with ≥1% prevalence
#   - Handles both MGBB-LLF and MGBB-ABC cohorts

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
})

# ---------- Command-line arguments ----------
option_list <- list(
  make_option("--cohort", type="character", help="Cohort name: MGBB-LLF or MGBB-ABC"),
  make_option("--virsight_promax", type="character", help="VirSIGHT Promax Hits CSV"),
  make_option("--husight_fl", type="character", help="HuSIGHT Full-Length Hits CSV"),
  make_option("--out_dir", type="character", default="results/binaries",
              help="Output directory [default %default]"),
  make_option("--min_prev_pct", type="double", default=1,
              help="Minimum prevalence (%) to retain a feature [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- Configuration ----------
if (opt$cohort == "MGBB-LLF") {
  sample_prefix <- "AM_"
  cohort_n <- 1290
} else if (opt$cohort == "MGBB-ABC") {
  sample_prefix <- "^10"
  cohort_n <- 300
} else {
  stop("Unknown cohort: must be MGBB-LLF or MGBB-ABC")
}

# HSV species reference
HSV_SPECIES <- c(
  "Human alphaherpesvirus 1", "Human alphaherpesvirus 2", "Human alphaherpesvirus 3",
  "Human betaherpesvirus 5", "Human betaherpesvirus 6A", "Human betaherpesvirus 6B",
  "Human betaherpesvirus 7", "Human gammaherpesvirus 4", "Human gammaherpesvirus 8"
)

# ---------- Helper functions ----------
bin_from_fob <- function(x) as.integer(!is.na(as.numeric(x)) & as.numeric(x) > 1)

apply_min_prevalence <- function(df, min_prev_pct) {
  n <- nrow(df)
  min_sum <- ceiling((min_prev_pct / 100) * n)
  keep <- vapply(df, function(x) sum(x, na.rm = TRUE) >= min_sum, logical(1))
  df[, keep, drop = FALSE]
}

# ---------- Process VirSIGHT HSV ----------
message("Reading VirSIGHT Promax file: ", opt$virsight_promax)
vi <- fread(opt$virsight_promax)
vi_hsv <- vi[taxon_species %in% HSV_SPECIES]

if (nrow(vi_hsv) == 0) stop("No HSV rows found in VirSIGHT data.")

sample_cols_vi <- grep(sample_prefix, names(vi_hsv), value = TRUE)

if (length(sample_cols_vi) == 0) stop("No usable virus sample columns found after removing controls.")

# Binary transform (>1 → 1)
mat_vi <- as.data.frame(vi_hsv[, ..sample_cols_vi])
mat_vi[] <- lapply(mat_vi, bin_from_fob)
mat_vi1 <- as.data.frame(t(mat_vi))
colnames(mat_vi1) <- vi_hsv$UniProt_acc
rownames(mat_vi1) <- NULL
mat_vi1 <- apply_min_prevalence(mat_vi1, opt$min_prev_pct)
mat_vi1$Sample_ID <- colnames(mat_vi)

# ---------- Process HuSIGHT Full-Length ----------
message("Reading HuSIGHT Full-Length file: ", opt$husight_fl)
hu <- fread(opt$husight_fl)
if (!"var_id" %in% names(hu)) hu[, var_id := paste0("h", .I)]

sample_cols_hu <- grep(sample_prefix, names(hu), value = TRUE)
if (length(sample_cols_hu) == 0) stop("No usable human sample columns found after removing controls.")

mat_hu <- as.data.frame(hu[, ..sample_cols_hu])
mat_hu[] <- lapply(mat_hu, bin_from_fob)
mat_hu1 <- as.data.frame(t(mat_hu))
colnames(mat_hu1) <- hu$var_id
rownames(mat_hu1) <- NULL
mat_hu1 <- apply_min_prevalence(mat_hu1, opt$min_prev_pct)
mat_hu1$Sample_ID <- colnames(mat_hu)

# ---------- Write Outputs ----------
out_virus <- file.path(opt$out_dir, sprintf("%s_hsv_binary.tsv", opt$cohort))
out_human <- file.path(opt$out_dir, sprintf("%s_human_fl_binary.tsv", opt$cohort))
fwrite(as.data.table(mat_vi1), out_virus, sep = "\t")
fwrite(as.data.table(mat_hu1), out_human, sep = "\t")

message("✔ Done building binary matrices for ", opt$cohort)
message("  - Virus HSV binary: ", normalizePath(out_virus))
message("  - Human FL binary:  ", normalizePath(out_human))
