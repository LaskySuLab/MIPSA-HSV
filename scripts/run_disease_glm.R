#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(dplyr); library(stringr); library(purrr); library(tidyr)
  source("R/utils_common.R"); source("R/utils_qtl.R")
})

# ----------------------------- CLI -----------------------------
option_list <- list(
  make_option("--phe",      type="character", help="MIPSA_Asthma_1290.csv"),
  make_option("--vi_bin",   type="character", help="hsv_proteins_binary.txt"),
  make_option("--hu_bin",   type="character", help="human_fl_binary.txt"),
  make_option("--dx_count", type="character",
              help="a directory containing prevalence/incidence files"),
  make_option("--out_dir",  type="character", default="results/Res_disease", help="Output dir")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------- IO ------------------------------
# Phenotype table and covariate processing
phe <- fread(opt$pheno)
phe[, Subject_Id := as.character(Subject_Id)]
phe <- build_disease_combine(phe)
phe <- recode_covariates(phe)

# Binary matrices
vi_bin <- fread(opt$vi_bin)  # virus peptides: columns = Sample_ID + UniProt_acc
hu_bin <- fread(opt$hu_bin)  # human FL:      columns = Sample_ID + hX

# Merge (right_join keeps only subjects present in the binary matrices)
merge_dt <- phe %>%
  right_join(hu_bin, by = "Sample_ID") %>%
  right_join(vi_bin, by = "Sample_ID") %>%
  as.data.table()

# ----------------------------- Disease cols -----------------------------
dx_cols <- grep("_combine$", colnames(merge_dt), value = TRUE)
if (length(dx_cols) == 0) stop("No *_combine disease columns found. Did build_disease_combine() run correctly?")

# Keep diseases with at least 13 positives
dx_sums <- colSums(merge_dt[, ..dx_cols], na.rm = TRUE)
dx_list <- names(dx_sums[dx_sums >= 13])

# Coerce disease columns to binary 0/1 integer (avoid factors)
merge_dt[, (dx_list) := lapply(.SD, function(z) as.integer(z != 0 & !is.na(z))), .SDcols = dx_list]

# ----------------------------- Run GLMs -----------------------------
# Human FL predictors (all antibodies)
rep_llf = fread(opt$rep_llf)
hu_genes <- unique(rep_llf$antibody)

# Function for parallel processing for human antibodies
run_glm_human <- function(ii) {
  dx <- dx_list[ii]
  print(dx)
  glm_results <- lapply(hu_genes, function(var) {
    formula <- as.formula(paste0(dx, " ~ ", var, " + Age_at_collect + Gender + BMI_most_recent + Race_Group + Smoking + Alcohol"))
    res <- tryCatch(glm(formula, family = "binomial", data = merge_dt), error = function(e) NULL)
    if (!is.null(res)) {
      summary_res <- summary(res)
      if (length(summary_res$coefficients) == 32) {
        return(c(summary_res$coefficients[2, 1:4]))
      }
    }
    return(rep(NA, 4))
  })
  outcome <- do.call(rbind, glm_results)
  colnames(outcome) <- c("Beta", "SE", "Z", "P")
  outcome <- as.data.table(outcome)
  outcome[, var_id := hu_genes]
  outcome[, disease := dx]
  
  fwrite(outcome[order(P)], file = file.path(opt$out_dir, paste0("bin_fchange_hu_",dx, "_com_glm.csv")), row.names = FALSE, col.names = TRUE, quote = FALSE)
}

lapply(1:length(dx_list), run_glm_human)


# Viral peptide predictors (all viral UniProt_acc)
virus_p <- setdiff(colnames(vi_bin), "Sample_ID")

# Function for parallel processing
run_glm_virus <- function(ii) {
  dx <- dx_list[ii]
  print(dx)
  glm_results <- lapply(virus_p, function(var) {
    formula <- as.formula(paste0(dx, " ~ ", var, " + Age_at_collect + Gender + BMI_most_recent + Race_Group + Smoking + Alcohol"))
    res <- tryCatch(glm(formula, family = "binomial", data = merge_dt), error = function(e) NULL)
    if (!is.null(res)) {
      summary_res <- summary(res)
      if (length(summary_res$coefficients) == 32) {
        return(c(summary_res$coefficients[2, 1:4]))
      }
    }
    return(rep(NA, 4))
  })
  outcome <- do.call(rbind, glm_results)
  colnames(outcome) <- c("Beta", "SE", "Z", "P")
  outcome <- as.data.table(outcome)
  outcome[, var_id := virus_p]
  outcome[, disease := dx]
  
  fwrite(outcome[order(P)], file = file.path(opt$out_dir,paste0("bin_fchange_vi_",dx, "_com_glm.csv")), row.names = FALSE, col.names = TRUE, quote = FALSE)
}

lapply(1:length(dx_list), run_glm_virus)


# ----------------------------- Combine per-DX CSVs -----------------------------
bind_csvs <- function(pattern, out_path) {
  files <- list.files(opt$out_dir, pattern = pattern, full.names = TRUE)
  if (!length(files)) return(invisible(NULL))
  DTs <- lapply(files, fread)
  all <- rbindlist(DTs, use.names = TRUE, fill = TRUE)
  fwrite(all, out_path)
  invisible(all)
}

human_all_path <- file.path(opt$out_dir, "tot_bin_fchange_FL_dx_com_glm.csv")
virus_all_path <- file.path(opt$out_dir, "tot_bin_fchange_Virus_dx_com_glm.csv")

human.res <- bind_csvs("^bin_fchange_hu_.*_com_glm\\.csv$", human_all_path)
virus.res <- bind_csvs("^bin_fchange_vi_.*_com_glm\\.csv$", virus_all_path)

if (is.null(human.res) || is.null(virus.res)) {
  message("Note: some combined outputs are empty—check input matrices and filters.")
}

# ----------------------------- Human results annotate & write -----------------------------
hu_fl <- fread(opt$hu_fl)
if (!"var_id" %in% names(hu_fl)) hu_fl[, var_id := paste0("h", seq_len(.N))]

if (!is.null(human.res) && nrow(human.res)) {
  human.res1 <- human.res %>%
    mutate(direction_of_effect = ifelse(Beta > 0, "positive", "negative"),
           disease = gsub("_combine$", "", disease)) %>%
    left_join(hu_fl[, .(gene_symbol, UniProt_acc, var_id)], by = c("var_id"))

    fwrite(human.res1,
         file.path(opt$out_dir, "rep_bin_fchange_FL_dx_com_glm_ann.tsv"),
         sep = "\t", quote = FALSE)
}

# ----------------------------- Virus results annotate & write -----------------------------
vi_promax <- fread(opt$vi_promax)

if (!is.null(virus.res) && nrow(virus.res)) {
  virus.res1 <- virus.res %>%
    mutate(direction_of_effect = ifelse(Beta > 0, "positive", "negative"),
           disease = gsub("_combine$", "", disease)) %>%
    left_join(vi_promax[, .(UniProt_acc, taxon_genus, taxon_species, product)],
              by = c("var_id" = "UniProt_acc"))

  fwrite(virus.res1,
         file.path(opt$out_dir, "rep_bin_fchange_HSV_dx_com_glm_ann.tsv"),
         sep = "\t", quote = FALSE)
}

message("Done. Outputs in: ", opt$out_dir)
