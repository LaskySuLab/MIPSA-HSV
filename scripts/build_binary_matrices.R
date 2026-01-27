#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  source("R/utils_common.R"); source("R/utils_qtl.R")
})

# ---------- Command-line arguments ----------
option_list <- list(
  make_option("--cohort", type="character", help="Cohort name: MGBB-LLF or MGBB-LEC"),
  make_option("--virsight_promax", type="character", help="VirSIGHT Promax Hits"),
  make_option("--husight_fl", type="character", help="HuSIGHT Full-Length Hits"),
  make_option("--out_dir", type="character", default="./Data",
              help="Output directory [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- Configuration ----------
# HSV species reference
HSV_SPECIES <- c(
  "Human alphaherpesvirus 1", "Human alphaherpesvirus 2", "Human alphaherpesvirus 3",
  "Human betaherpesvirus 5", "Human betaherpesvirus 6A", "Human betaherpesvirus 6B",
  "Human betaherpesvirus 7", "Human gammaherpesvirus 4", "Human gammaherpesvirus 8"
)

# ---------- Helper functions ----------
bin_from_fob <- function(x) as.integer(!is.na(as.numeric(x)) & as.numeric(x) > 1)

apply_min_prevalence <- function(df, min_sum) {
  n <- nrow(df)
  keep <- vapply(df, function(x) sum(x, na.rm = TRUE) > min_sum, logical(1))
  df[, keep, drop = FALSE]
}

# ----------LLF: Process VirSIGHT HSV ----------
pheno.llf <- fread(opt$pheno.llf)

vi.llf <- fread(opt$virsight_promax.llf)
sample_lists = pheno.llf$Sample_ID

vi_hsv <- vi.llf[taxon_species %in% HSV_SPECIES]
if (nrow(vi_hsv) == 0) stop("No HSV rows found in VirSIGHT data.")

# Binary transform (>1 → 1)
mat_vi <- as.data.frame(subset(vi_hsv, select = sample_lists))
mat_vi[] <- lapply(mat_vi, bin_from_fob)
mat_vi1 <- as.data.frame(t(mat_vi))
colnames(mat_vi1) <- vi_hsv$UniProt_acc
rownames(mat_vi1) <- NULL
mat_vi1 <- apply_min_prevalence(mat_vi1, length(sample_lists)/100)
dim(mat_vi1)
# [1] 1290 2323
mat_vi1$Sample_ID <- colnames(mat_vi)

# ----------LLF: Process HuSIGHT Full-Length ----------
hu <- fread(opt$husight_fl.llf)
if (!"var_id" %in% names(hu)) hu[, var_id := paste0("h", .I)]

mat_hu <- as.data.frame(subset(hu, select = sample_lists))
mat_hu[] <- lapply(mat_hu, bin_from_fob)
mat_hu1 <- as.data.frame(t(mat_hu))
colnames(mat_hu1) <- hu$var_id
rownames(mat_hu1) <- NULL
mat_hu1 <- apply_min_prevalence(mat_hu1, length(sample_lists)/100)
dim(mat_hu1)
# [1] 1290 2258
mat_hu1$Sample_ID <- colnames(mat_hu)

# ----------LLF: Write Outputs ----------
out_virus <- file.path(opt$out_dir, sprintf("hsv_promax_bin_%s.tsv", 'MGBB-LLF'))
out_human <- file.path(opt$out_dir, sprintf("human_fl_bin_%s.tsv", 'MGBB-LLF'))
fwrite(as.data.table(mat_vi1), out_virus, sep = "\t")
fwrite(as.data.table(mat_hu1), out_human, sep = "\t")




# ----------LEC: Process VirSIGHT HSV ----------
pheno.lec <- fread(opt$pheno.lec)
sample_lists = pheno.lec$Sample_Id

vi.abc <- fread(opt$virsight_promax.abc)
vi.leo <- fread(opt$virsight_promax.leo)
all(vi.abc$UniProt_acc==vi.leo$UniProt_acc)
# [1] TRUE

vi.lec = cbind(vi.abc, vi.leo)

vi_hsv <- vi.lec[taxon_species %in% HSV_SPECIES]

# Binary transform (>1 → 1)
mat_vi <- as.data.frame(subset(vi_hsv, select = sample_lists))
mat_vi[] <- lapply(mat_vi, bin_from_fob)
mat_vi1 <- as.data.frame(t(mat_vi))
colnames(mat_vi1) <- vi_hsv$UniProt_acc
rownames(mat_vi1) <- NULL
mat_vi1 <- apply_min_prevalence(mat_vi1, length(sample_lists)/100)
dim(mat_vi1)
# [1]  763 2342
mat_vi1$Sample_ID <- colnames(mat_vi)

# ----------LEC: Process HuSIGHT Full-Length ----------
hu.abc <- fread(opt$husight_fl.abc)
if (!"var_id" %in% names(hu.abc)) hu.abc[, var_id := paste0("h", .I)]

hu.leo <- fread(opt$husight_fl.leo)
hu.abc1 = subset(hu.abc, hu.abc$pep_id%in%hu.leo$pep_id)
all(hu.leo$pep_id==hu.abc1$pep_id)

hu.lec = cbind(hu.abc1, hu.leo)

mat_hu <- as.data.frame(subset(hu.lec, select = sample_lists))
mat_hu[] <- lapply(mat_hu, bin_from_fob)
mat_hu1 <- as.data.frame(t(mat_hu))
colnames(mat_hu1) <- hu.abc1$var_id
rownames(mat_hu1) <- NULL
mat_hu1 <- apply_min_prevalence(mat_hu1, length(sample_lists)/100)
dim(mat_hu1)
# [1]  763 1064
mat_hu1$Sample_ID <- colnames(mat_hu)

# ----------LEC: Write Outputs ----------
out_virus <- file.path(opt$out_dir, sprintf("hsv_promax_bin_%s.tsv", 'MGBB-LEC'))
out_human <- file.path(opt$out_dir, sprintf("human_fl_bin_%s.tsv", 'MGBB-LEC'))
fwrite(as.data.table(mat_vi1), out_virus, sep = "\t")
fwrite(as.data.table(mat_hu1), out_human, sep = "\t")

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  source("R/utils_common.R"); source("R/utils_qtl.R")
})

# ---------- Command-line arguments ----------
option_list <- list(
  make_option("--cohort", type="character", help="Cohort name: MGBB-LLF or MGBB-LEC"),
  make_option("--virsight_promax", type="character", help="VirSIGHT Promax Hits"),
  make_option("--husight_fl", type="character", help="HuSIGHT Full-Length Hits"),
  make_option("--out_dir", type="character", default="./Data",
              help="Output directory [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- Configuration ----------
# HSV species reference
HSV_SPECIES <- c(
  "Human alphaherpesvirus 1", "Human alphaherpesvirus 2", "Human alphaherpesvirus 3",
  "Human betaherpesvirus 5", "Human betaherpesvirus 6A", "Human betaherpesvirus 6B",
  "Human betaherpesvirus 7", "Human gammaherpesvirus 4", "Human gammaherpesvirus 8"
)

# ---------- Helper functions ----------
bin_from_fob <- function(x) as.integer(!is.na(as.numeric(x)) & as.numeric(x) > 1)

apply_min_prevalence <- function(df, min_sum) {
  n <- nrow(df)
  keep <- vapply(df, function(x) sum(x, na.rm = TRUE) > min_sum, logical(1))
  df[, keep, drop = FALSE]
}

# ----------LLF: Process VirSIGHT HSV ----------
pheno.llf <- fread(opt$pheno.llf)

vi.llf <- fread(opt$virsight_promax.llf)
sample_lists = pheno.llf$Sample_ID

vi_hsv <- vi.llf[taxon_species %in% HSV_SPECIES]
if (nrow(vi_hsv) == 0) stop("No HSV rows found in VirSIGHT data.")

# Binary transform (>1 → 1)
mat_vi <- as.data.frame(subset(vi_hsv, select = sample_lists))
mat_vi[] <- lapply(mat_vi, bin_from_fob)
mat_vi1 <- as.data.frame(t(mat_vi))
colnames(mat_vi1) <- vi_hsv$UniProt_acc
rownames(mat_vi1) <- NULL
mat_vi1 <- apply_min_prevalence(mat_vi1, length(sample_lists)/100)
dim(mat_vi1)
# [1] 1290 2323
mat_vi1$Sample_ID <- colnames(mat_vi)

# ----------LLF: Process HuSIGHT Full-Length ----------
hu <- fread(opt$husight_fl.llf)
if (!"var_id" %in% names(hu)) hu[, var_id := paste0("h", .I)]

mat_hu <- as.data.frame(subset(hu, select = sample_lists))
mat_hu[] <- lapply(mat_hu, bin_from_fob)
mat_hu1 <- as.data.frame(t(mat_hu))
colnames(mat_hu1) <- hu$var_id
rownames(mat_hu1) <- NULL
mat_hu1 <- apply_min_prevalence(mat_hu1, length(sample_lists)/100)
dim(mat_hu1)
# [1] 1290 2258
mat_hu1$Sample_ID <- colnames(mat_hu)

# ----------LLF: Write Outputs ----------
out_virus <- file.path(opt$out_dir, sprintf("hsv_promax_bin_%s.tsv", 'MGBB-LLF'))
out_human <- file.path(opt$out_dir, sprintf("human_fl_bin_%s.tsv", 'MGBB-LLF'))
fwrite(as.data.table(mat_vi1), out_virus, sep = "\t")
fwrite(as.data.table(mat_hu1), out_human, sep = "\t")




# ----------LEC: Process VirSIGHT HSV ----------
pheno.lec <- fread(opt$pheno.lec)
sample_lists = pheno.lec$Sample_Id

vi.abc <- fread(opt$virsight_promax.abc)
vi.leo <- fread(opt$virsight_promax.leo)
all(vi.abc$UniProt_acc==vi.leo$UniProt_acc)
# [1] TRUE

vi.lec = cbind(vi.abc, vi.leo)

vi_hsv <- vi.lec[taxon_species %in% HSV_SPECIES]

# Binary transform (>1 → 1)
mat_vi <- as.data.frame(subset(vi_hsv, select = sample_lists))
mat_vi[] <- lapply(mat_vi, bin_from_fob)
mat_vi1 <- as.data.frame(t(mat_vi))
colnames(mat_vi1) <- vi_hsv$UniProt_acc
rownames(mat_vi1) <- NULL
mat_vi1 <- apply_min_prevalence(mat_vi1, length(sample_lists)/100)
dim(mat_vi1)
# [1]  763 2342
mat_vi1$Sample_ID <- colnames(mat_vi)

# ----------LEC: Process HuSIGHT Full-Length ----------
hu.abc <- fread(opt$husight_fl.abc)
if (!"var_id" %in% names(hu.abc)) hu.abc[, var_id := paste0("h", .I)]

hu.leo <- fread(opt$husight_fl.leo)
hu.abc1 = subset(hu.abc, hu.abc$pep_id%in%hu.leo$pep_id)
all(hu.leo$pep_id==hu.abc1$pep_id)

hu.lec = cbind(hu.abc1, hu.leo)

mat_hu <- as.data.frame(subset(hu.lec, select = sample_lists))
mat_hu[] <- lapply(mat_hu, bin_from_fob)
mat_hu1 <- as.data.frame(t(mat_hu))
colnames(mat_hu1) <- hu.abc1$var_id
rownames(mat_hu1) <- NULL
mat_hu1 <- apply_min_prevalence(mat_hu1, length(sample_lists)/100)
dim(mat_hu1)
# [1]  763 1064
mat_hu1$Sample_ID <- colnames(mat_hu)

# ----------LEC: Write Outputs ----------
out_virus <- file.path(opt$out_dir, sprintf("hsv_promax_bin_%s.tsv", 'MGBB-LEC'))
out_human <- file.path(opt$out_dir, sprintf("human_fl_bin_%s.tsv", 'MGBB-LEC'))
fwrite(as.data.table(mat_vi1), out_virus, sep = "\t")
fwrite(as.data.table(mat_hu1), out_human, sep = "\t")

