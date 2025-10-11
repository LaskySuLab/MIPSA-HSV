#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(dplyr)
  source("R/utils_common.R")
})

opt <- list()
option_list <- list(
  make_option("--vi_csv",  type="character", help="VirSIGHT Promax FoB CSV"),
  make_option("--hu_csv",  type="character", help="HuSIGHT FullLength FoB CSV"),
  make_option("--out_dir", type="character", default="results", help="Output dir")
)
opt <- parse_args(OptionParser(option_list=option_list))

ensure_dir(opt$out_dir)

# Read raw FoB and convert to binary 0/1 (non-zero -> 1)
vi <- fread(opt$vi_csv)
hu <- fread(opt$hu_csv)

# Virus: wide binary with Subject_Id rows
vi_w <- as.data.frame(t(vi[, 6:ncol(vi)]))
colnames(vi_w) <- vi$UniProt_acc
vi_w[vi_w == 1] <- 0
vi_w[vi_w != 0] <- 1
vi_w$Subject_Id <- rownames(vi_w)

# Human: assign var_id (h1..hN) if absent; wide binary
if (!"var_id" %in% names(hu)) hu$var_id <- paste0("h", seq_len(nrow(hu)))
hu_w <- as.data.frame(t(hu[, 11:ncol(hu)]))
colnames(hu_w) <- hu$var_id
hu_w[hu_w == 1] <- 0
hu_w[hu_w != 0] <- 1
hu_w$Subject_Id <- rownames(hu_w)

# Save TSVs
fwrite(as.data.table(vi_w), file.path(opt$out_dir, "virus_proteins_binary.txt"), sep="\t")
fwrite(as.data.table(hu_w), file.path(opt$out_dir, "human_fl_binary.txt"), sep="\t")

# Optional log (if needed later)
vi_log <- as.data.frame(t(vi[, 6:ncol(vi)])); colnames(vi_log) <- vi$UniProt_acc
hu_log <- as.data.frame(t(hu[, 11:ncol(hu)])); colnames(hu_log) <- hu$var_id
vi_log$Subject_Id <- rownames(vi_log); hu_log$Subject_Id <- rownames(hu_log)
fwrite(as.data.table(vi_log), file.path(opt$out_dir, "hsv_proteins_log.txt"), sep="\t")
fwrite(as.data.table(hu_log), file.path(opt$out_dir, "human_fl_log.txt"), sep="\t")
