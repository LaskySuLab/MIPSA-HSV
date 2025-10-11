#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(dplyr); library(stringr); library(purrr)
  source("R/utils_common.R"); source("R/utils_qtl.R")
})

option_list <- list(
  make_option("--phe",     type="character", help="MIPSA_Asthma_1290.csv"),
  make_option("--vi_bin",  type="character", help="hsv_proteins_binary.txt"),
  make_option("--hu_bin",  type="character", help="human_fl_binary.txt"),
  make_option("--dx_prev", type="character", help="dx_prevalence_sums"),
  make_option("--dx_inc",  type="character", help="dx_incidence_sums"),
  make_option("--out_dir", type="character", default="results/Res_disease", help="Output dir")
)
opt <- parse_args(OptionParser(option_list=option_list))
ensure_dir(opt$out_dir)

phe <- fread(opt$phe); phe[, Subject_Id := as.character(Subject_Id)]
phe <- build_disease_combine(phe)
phe <- recode_covariates(phe)

vi_bin <- fread(opt$vi_bin)
hu_bin <- fread(opt$hu_bin)

merge_dt <- phe %>%
  right_join(hu_bin, by=c("Sample_ID"="Subject_Id")) %>%
  right_join(vi_bin, by=c("Sample_ID"="Subject_Id")) %>% as.data.table()

# disease columns
dx_cols <- grep("_combine$", colnames(merge_dt), value=TRUE)
dx_status <- merge_dt[, ..dx_cols]
dx_status <- dx_status[, which(colSums(.SD, na.rm=TRUE) >= 13), .SDcols = dx_cols, with = FALSE]
dx_list <- colnames(dx_status)

covar <- c("Age_at_collect","Gender","BMI_most_recent","Race_Group","Smoking","Alcohol")
stopif_missing_cols(merge_dt, covar)

# Hu predictors
hu_genes <- colnames(hu_bin)[colnames(hu_bin)!="Subject_Id"]
for (dx in dx_list) {
  res <- map_dfr(hu_genes, function(hg) {
    if (length(unique(merge_dt[[hg]])) < 2) return(NULL)
    fml <- as.formula(paste0(dx, " ~ ", hg, " + ", paste(covar, collapse=" + ")))
    s <- glm_safe(fml, merge_dt)
    if (is.null(s)) return(NULL)
    cbind(data.frame(disease=dx, var_id=hg, type="human"), s)
  })
  if (nrow(res)) fwrite(res[order(P)], file.path(opt$out_dir, paste0("bin_fchange_hu_", dx, "_com_glm.csv")))
}

# Virus predictors
virus_p <- setdiff(colnames(vi_bin), "Subject_Id")
for (dx in dx_list) {
  res <- map_dfr(virus_p, function(vp) {
    if (length(unique(merge_dt[[vp]])) < 2) return(NULL)
    fml <- as.formula(paste0(dx, " ~ ", vp, " + ", paste(covar, collapse=" + ")))
    s <- glm_safe(fml, merge_dt)
    if (is.null(s)) return(NULL)
    cbind(data.frame(disease=dx, var_id=vp, type="virus"), s)
  })
  if (nrow(res)) fwrite(res[order(P)], file.path(opt$out_dir, paste0("bin_fchange_vi_", dx, "_com_glm.csv")))
}
