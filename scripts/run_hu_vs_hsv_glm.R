#!/usr/bin/env Rscript
# Run Hu (FL) ~ HSV (Promax) GLMs for LLF and ABC cohorts, annotate, FDR, and replicate
# Outputs:
#   - results/hsv_qtl/llf_hsv_bin_glm_all.tsv
#   - results/hsv_qtl/llf_hsv_bin_glm_sig.tsv
#   - results/hsv_qtl/abc_hsv_bin_glm_all.tsv
#   - results/hsv_qtl/abc_hsv_bin_glm_sig.tsv
#   - results/hsv_qtl/rep_hsv.sig              (LLF significant replicated in ABC)
#   - results/hsv_qtl/abc_hsv.rep             (matching ABC rows for the replicated set)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(parallel)
})

# ---------------------------- CLI ----------------------------
opt_list <- list(
  # LLF inputs
  make_option("--llf_phe", type="character",
              help="LLF phenotype CSV (e.g., MIPSA_Asthma_1290.csv)"),
  make_option("--llf_vi_promax", type="character",
              help="LLF VirSIGHT Promax CSV for annotation"),
  make_option("--llf_hu_fl", type="character",
              help="LLF HuSIGHT FullLength CSV for annotation (h1..hN)"),
  make_option("--llf_vi_bin", type="character",
              help="LLF HSV binary TSV (from build_binary_matrices.R)"),
  make_option("--llf_hu_bin", type="character",
              help="LLF human FL binary TSV (from build_binary_matrices.R)"),
  
  # ABC inputs
  make_option("--abc_phe", type="character",
              help="ABC phenotype CSV (e.g., Ab_pheno.csv)"),
  make_option("--abc_vi_promax", type="character",
              help="ABC VirSIGHT Promax CSV for annotation"),
  make_option("--abc_hu_fl", type="character",
              help="ABC HuSIGHT FullLength CSV for annotation (h1..hN)"),
  make_option("--abc_vi_bin", type="character",
              help="ABC HSV binary TSV (from build_binary_matrices.R)"),
  make_option("--abc_hu_bin", type="character",
              help="ABC human FL binary TSV (from build_binary_matrices.R)"),
  
  # General
  make_option("--out_dir", type="character", default="results/hsv_qtl",
              help="Output directory [default %default]"),
  make_option("--n_cores", type="integer", default=max(1, detectCores()-1),
              help="CPU cores to use [default %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------- Helpers --------------------------
recode_covariates <- function(dt) {
  # Race_Group: White -> 0, others -> 1
  if ("Race_Group" %in% names(dt)) {
    dt$Race_Group <- ifelse(dt$Race_Group == "White", 0, 1)
  }
  # Gender: M -> 0, others -> 1
  if ("Gender" %in% names(dt)) {
    dt$Gender <- ifelse(dt$Gender == "M", 0, 1)
  }
  
  # Smoking (combine biobank + RPDR; unknown fallback)
  smoking <- c("Current-smoker", "Former-smoker")
  if ("Smoking_flag_Biobank" %in% names(dt)) {
    dt[, Smoking_flag_Biobank := ifelse(Smoking_flag_Biobank %in% smoking, "Current/Former-smoker", Smoking_flag_Biobank)]
    dt[Smoking_flag_Biobank == "", Smoking_flag_Biobank := NA]
  }
  if ("Smoking_RPDR_most_recent" %in% names(dt)) {
    dt[, Smoking_RPDR_most_recent := ifelse(Smoking_RPDR_most_recent %in% smoking, "Current/Former-smoker", Smoking_RPDR_most_recent)]
    dt[Smoking_RPDR_most_recent == "", Smoking_RPDR_most_recent := NA]
  }
  dt[, Smoking_flag := if ("Smoking_flag_Biobank" %in% names(dt)) Smoking_flag_Biobank else NA_character_]
  if ("Smoking_RPDR_most_recent" %in% names(dt)) {
    dt[is.na(Smoking_flag), Smoking_flag := Smoking_RPDR_most_recent]
  }
  dt[is.na(Smoking_flag), Smoking_flag := "Unknown"]
  dt[, Smoking := ifelse(Smoking_flag %in% c("Current/Former-smoker"), 1, 0)]
  
  # Alcohol
  if ("Alcohol_flag_Biobank" %in% names(dt)) {
    alc_levels <- sort(unique(dt$Alcohol_flag_Biobank))
    alc_levels <- alc_levels[!is.na(alc_levels)]
    if (length(alc_levels) >= 2) {
      # cohort-specific slices will differ; we apply safe recoding:
      dt[, Alcohol_flag_Biobank := ifelse(Alcohol_flag_Biobank %in% c("Current-drinker", "Former-drinker"), "Current/Former-drinker", Alcohol_flag_Biobank)]
      dt[, Alcohol_flag_Biobank := ifelse(Alcohol_flag_Biobank %in% c(" 1. None, or less than 1 per month "), "No-alcohol", Alcohol_flag_Biobank)]
      dt[Alcohol_flag_Biobank == "", Alcohol_flag_Biobank := NA]
    }
  }
  if ("Alcohol_RPDR_most_recent" %in% names(dt)) {
    dt[, Alcohol_RPDR_most_recent := ifelse(Alcohol_RPDR_most_recent %in% c("Current-drinker", "Former-drinker"), "Current/Former-drinker", Alcohol_RPDR_most_recent)]
    dt[Alcohol_RPDR_most_recent == "", Alcohol_RPDR_most_recent := NA]
  }
  dt[, Alcohol_flag := if ("Alcohol_flag_Biobank" %in% names(dt)) Alcohol_flag_Biobank else NA_character_]
  if ("Alcohol_RPDR_most_recent" %in% names(dt)) {
    dt[is.na(Alcohol_flag), Alcohol_flag := Alcohol_RPDR_most_recent]
  }
  dt[is.na(Alcohol_flag), Alcohol_flag := "Unknown"]
  dt[, Alcohol := ifelse(Alcohol_flag %in% c("Current/Former-drinker"), 1, 0)]
  
  dt
}

glm_block <- function(df, hu_var, vi_var, covars) {
  # Returns one data.table with columns: Beta, SE, Z, P, var_id, antibody
  # antibody ~ vi + covars
  glm_results <- lapply(vi_var, function(vv) {
    ftxt <- paste(hu_var, "~", paste(c(vv, covars), collapse=" + "))
    res <- tryCatch(glm(as.formula(ftxt), family = "binomial", data = df), error = function(e) NULL)
    if (is.null(res)) return(c(NA, NA, NA, NA))
    sm <- summary(res)
    if (nrow(sm$coefficients) < 2) return(c(NA, NA, NA, NA))
    sm$coefficients[2, 1:4]
  })
  out <- as.data.table(do.call(rbind, glm_results))
  setnames(out, c("Beta","SE","Z","P"))
  out[, var_id := vi_var]
  out[, antibody := hu_var]
  out[]
}

run_glm_parallel <- function(df, hu_vars, vi_vars, covars, n_cores) {
  n_cores <- max(1, n_cores)
  if (n_cores == 1) {
    outs <- lapply(hu_vars, function(hh) glm_block(df, hh, vi_vars, covars))
  } else {
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl))
    clusterExport(cl, varlist = c("df","hu_vars","vi_vars","covars","glm_block"),
                  envir = environment())
    clusterEvalQ(cl, { library(data.table); library(stats) })
    outs <- parLapply(cl, hu_vars, function(hh) glm_block(df, hh, vi_vars, covars))
  }
  rbindlist(outs, use.names = TRUE, fill = TRUE)
}

annotate_llf <- function(tab, vi_promax, hu_fl) {
  # LLF annotations
  hu_fl$antibody <- paste0("h", seq_len(nrow(hu_fl)))
  tab %>%
    left_join(hu_fl[, c("gene_symbol","UniProt_acc","antibody")], by = c("antibody")) %>%
    left_join(vi_promax[, c("UniProt_acc","taxon_genus","taxon_species","product")],
              by = c("var_id" = "UniProt_acc")) %>%
    as.data.table()
}

annotate_abc <- function(tab, vi_promax, hu_fl) {
  hu_fl$var_id <- paste0("h", seq_len(nrow(hu_fl)))
  tab %>%
    left_join(vi_promax[, 1:5], by = c("var_id" = "UniProt_acc")) %>%
    left_join(hu_fl[, c("UniProt_acc","gene_symbol","var_id")], by = c("antibody" = "var_id")) %>%
    as.data.table()
}

# Key used for replication check
make_qtl_key <- function(dt) {
  dt$direction_of_effect <- ifelse(dt$Beta > 0, "positive", "negative")
  dt$qtl <- paste0(dt$var_id, "_", dt$antibody, "_", dt$direction_of_effect)
  dt
}

# -------------------------- LLF block --------------------------
if (!all(c(opt$llf_phe, opt$llf_vi_bin, opt$llf_hu_bin, opt$llf_vi_promax, opt$llf_hu_fl) %in% "") &&
    !any(is.null(c(opt$llf_phe, opt$llf_vi_bin, opt$llf_hu_bin, opt$llf_vi_promax, opt$llf_hu_fl)))) {
  
  message("=== LLF: reading inputs ===")
  llf_phe      <- fread(opt$llf_phe)
  llf_vi_bin   <- fread(opt$llf_vi_bin)     # columns: virus proteins + Sample_ID
  llf_hu_bin   <- fread(opt$llf_hu_bin)     # columns: human FL var_id + Sample_ID
  llf_vi_annot <- fread(opt$llf_vi_promax)
  llf_hu_annot <- fread(opt$llf_hu_fl)
  
  # Merge: phenotype + virus + human
  llf_merged <- llf_phe %>%
    right_join(llf_vi_bin, by = "Sample_ID") %>%
    right_join(llf_hu_bin, by = "Sample_ID") %>%
    as.data.table()
  
  llf_merged <- recode_covariates(llf_merged)
  
  # Variables
  hu_cols_llf <- setdiff(colnames(llf_hu_bin), "Sample_ID")
  vi_cols_llf <- setdiff(colnames(llf_vi_bin), "Sample_ID")
  covars_llf  <- c("Age_at_collect", "Gender", "BMI_most_recent", "Race_Group", "Smoking", "Alcohol")
  missing_cov <- setdiff(covars_llf, names(llf_merged))
  if (length(missing_cov)) stop("LLF: missing covariates: ", paste(missing_cov, collapse=", "))
  
  message("=== LLF: running GLMs in parallel ===")
  llf_res <- run_glm_parallel(
    df     = llf_merged,
    hu_vars= hu_cols_llf,
    vi_vars= vi_cols_llf,
    covars = covars_llf,
    n_cores= opt$n_cores
  )
  
  setorder(llf_res, P)
  # Annotate + FDR + direction
  llf_res <- annotate_llf(llf_res, vi_promax = llf_vi_annot, hu_fl = llf_hu_annot)
  llf_res$P.adj <- p.adjust(llf_res$P, method = "fdr")
  llf_res <- make_qtl_key(llf_res)
  
  # Save all + sig
  out_all_llf <- file.path(opt$out_dir, "llf_hsv_bin_glm_all.tsv")
  out_sig_llf <- file.path(opt$out_dir, "llf_hsv_bin_glm_sig.tsv")
  fwrite(llf_res, out_all_llf, sep = "\t")
  fwrite(llf_res[P.adj < 0.05], out_sig_llf, sep = "\t")
  
  # Report cutoff like your code
  sig_idx <- which(llf_res$P.adj <= 0.05)
  if (length(sig_idx)) {
    pval_cutoff <- max(llf_res$P[sig_idx], na.rm = TRUE)
    message(sprintf("LLF: FDR<=0.05 p-value cutoff ≈ %.10f", pval_cutoff))
  } else {
    message("LLF: No FDR-significant results at 0.05.")
  }
}

# -------------------------- ABC block --------------------------
if (!all(c(opt$abc_phe, opt$abc_vi_bin, opt$abc_hu_bin, opt$abc_vi_promax, opt$abc_hu_fl) %in% "") &&
    !any(is.null(c(opt$abc_phe, opt$abc_vi_bin, opt$abc_hu_bin, opt$abc_vi_promax, opt$abc_hu_fl)))) {
  
  message("=== ABC: reading inputs ===")
  abc_phe      <- fread(opt$abc_phe)
  abc_vi_bin   <- fread(opt$abc_vi_bin)     # columns: virus proteins + Sample_ID
  abc_hu_bin   <- fread(opt$abc_hu_bin)     # columns: human FL var_id + Sample_ID
  abc_vi_annot <- fread(opt$abc_vi_promax)
  abc_hu_annot <- fread(opt$abc_hu_fl)
  
  # Merge: phenotype uses Subject_Id; binaries use Sample_ID
  if (!("Subject_Id" %in% names(abc_phe))) {
    stop("ABC phenotype must contain 'Subject_Id' to join ABC binaries.")
  }
  abc_merged <- abc_phe %>%
    right_join(abc_vi_bin,  by = c("Subject_Id" = "Sample_ID")) %>%
    right_join(abc_hu_bin,  by = c("Subject_Id" = "Sample_ID")) %>%
    as.data.table()
  
  abc_merged <- recode_covariates(abc_merged)
  
  # Variables
  hu_cols_abc <- setdiff(colnames(abc_hu_bin), "Sample_ID")
  vi_cols_abc <- setdiff(colnames(abc_vi_bin), "Sample_ID")
  covars_abc  <- c("Age_at_collect", "Gender", "BMI_most_recent", "Race_Group", "Smoking", "Alcohol")
  missing_cov <- setdiff(covars_abc, names(abc_merged))
  if (length(missing_cov)) stop("ABC: missing covariates: ", paste(missing_cov, collapse=", "))
  
  message("=== ABC: running GLMs in parallel ===")
  abc_res <- run_glm_parallel(
    df      = abc_merged,
    hu_vars = hu_cols_abc,
    vi_vars = vi_cols_abc,
    covars  = covars_abc,
    n_cores = opt$n_cores
  )
  
  setorder(abc_res, P)
  # Annotate + FDR + direction
  abc_res <- annotate_abc(abc_res, vi_promax = abc_vi_annot, hu_fl = abc_hu_annot)
  abc_res$P.adj <- p.adjust(abc_res$P, method = "fdr")
  abc_res <- make_qtl_key(abc_res)
  
  # Save all + sig
  out_all_abc <- file.path(opt$out_dir, "abc_hsv_bin_glm_all.tsv")
  out_sig_abc <- file.path(opt$out_dir, "abc_hsv_bin_glm_sig.tsv")
  fwrite(abc_res, out_all_abc, sep = "\t")
  fwrite(abc_res[P.adj < 0.05], out_sig_abc, sep = "\t")
  
  sig_idx <- which(abc_res$P.adj <= 0.05)
  if (length(sig_idx)) {
    pval_cutoff <- max(abc_res$P[sig_idx], na.rm = TRUE)
    message(sprintf("ABC: FDR<=0.05 p-value cutoff ≈ %.10f", pval_cutoff))
  } else {
    message("ABC: No FDR-significant results at 0.05.")
  }
  
  # ------------------ Replication (LLF -> ABC) ------------------
  if (exists("llf_res")) {
    # LLF significant with gene symbols cleaned like your code
    llf_sig <- copy(llf_res[P.adj < 0.05])
    if ("gene_symbol" %in% names(llf_sig)) {
      llf_sig$gene_symbol <- ifelse(
        grepl("^TSPY3,TSPY10,LOC", llf_sig$gene_symbol),
        "TSPY3, TSPY10..",
        llf_sig$gene_symbol
      )
    }
    
    # ABC significant set
    abc_sig <- copy(abc_res[P.adj < 0.05])
    
    # Key already built: qtl = var_id_antibody_direction
    llf_keys <- unique(llf_sig$qtl)
    abc_keys <- unique(abc_sig$qtl)
    common   <- intersect(llf_keys, abc_keys)
    
    rep_llf <- llf_sig[qtl %in% common]
    rep_abc <- abc_sig[qtl %in% common]
    
    # Save exactly as your filenames
    fwrite(rep_llf, file.path(opt$out_dir, "rep_hsv.sig"), sep = "\t")
    fwrite(rep_abc, file.path(opt$out_dir, "abc_hsv.rep"), sep = "\t")
    message("Replication complete: N = ", length(common))
  } else {
    message("Replication skipped (LLF results not available in this run).")
  }
}

message("Done. Results written to: ", normalizePath(opt$out_dir))
