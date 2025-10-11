#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  source("R/utils_common.R")
  source("R/utils_qtl.R")
})

opt <- OptionParser(option_list = list(
  make_option("--cohort",         type="character", help="MGBB-LLF or MGBB-ABC"),
  make_option("--virsight-promax",type="character", help="VirSIGHT *Promax_Hits_Fold-Over-Background.csv"),
  make_option("--husight-fl",     type="character", help="HuSIGHT FullLength *Hits_Fold-Over-Background.csv"),
  make_option("--outdir",         type="character", default="results/binaries")
)) |> parse_args()

cfg <- cohort_preset(opt$cohort)
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

# ---- Virus HSV binary ----
vi <- fread(opt$virsight_promax)
vi_hsv <- vi[taxon_species %in% HSV_SPECIES]
# Keep only measurement columns
mat <- as.data.frame(vi_hsv[, starts_with(cfg$sample_prefix), with=FALSE])
mat[mat == 1] <- 0
mat[mat != 0] <- 1
mat <- as.data.frame(t(mat))
colnames(mat) <- vi_hsv$UniProt_acc

# 1% prevalence threshold
min_sum <- ceiling(0.01 * cfg$n)
mat_v <- mat %>% dplyr::select(where(~ is.numeric(.) && sum(.) > min_sum))
mat_v$Sample_ID <- rownames(mat_v)

fwrite(mat_v, file.path(opt$outdir, sprintf("%s_hsv_binary.tsv", opt$cohort)), sep="\t")

# ---- Human FL binary ----
hu <- fread(opt$husight_fl)
if (!"var_id" %in% names(hu)) hu[, var_id := paste0("h", .I)]
mat_h <- as.data.frame(hu[, starts_with(cfg$sample_prefix), with=FALSE])
mat_h[mat_h == 1] <- 0
mat_h[mat_h != 0] <- 1
mat_h <- as.data.frame(t(mat_h))
colnames(mat_h) <- hu$var_id

# 1% prevalence
mat_h2 <- mat_h %>% dplyr::select(where(~ is.numeric(.) && sum(.) > min_sum))
mat_h2$Sample_ID <- rownames(mat_h2)

fwrite(mat_h2, file.path(opt$outdir, sprintf("%s_human_fl_binary.tsv", opt$cohort)), sep="\t")

# LLF
Rscript scripts/build_binary_matrices.R \
  --cohort MGBB-LLF \
  --virsight-promax data/IB1007_VirSIGHT_Promax_Hits_Fold-Over-Background.csv \
  --husight-fl     data/IB1007_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv

# ABC
Rscript scripts/build_binary_matrices.R \
  --cohort MGBB-ABC \
  --virsight-promax data/IB1021_VirSIGHT_Promax_Hits_Fold-Over-Background.csv \
  --husight-fl     data/IB1021_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv

