#!/usr/bin/env Rscript
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
              help="llf_1290_phe.tsv"),
  make_option("--llf_vi_promax", type="character",
              help="IB1007_VirSIGHT_Promax_Hits_Fold-Over-Background.csv"),
  make_option("--llf_hu_fl", type="character",
              help="IB1007_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv"),
  make_option("--llf_vi_bin", type="character",
              help="hsv_promax_bin_MGBB-LLF.tsv"),
  make_option("--llf_hu_bin", type="character",
              help="human_fl_bin_MGBB-LLF.tsv"),
  
  # LEC inputs
  make_option("--lec_phe", type="character",
              help="lec_763_phe.tsv"),
  make_option("--lec_vi_bin", type="character",
              help="hsv_promax_bin_MGBB-LEC.tsv"),
  make_option("--lec_hu_bin", type="character",
              help="human_fl_bin_MGBB-LEC.tsv"),
  
  # General
  make_option("--out_dir", type="character", default="results/Run",
              help="Output directory [default %default]"),
  make_option("--n_cores", type="integer", default=max(1, detectCores()-1),
              help="CPU cores to use [default %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------- Helpers --------------------------
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

annotate_lec <- function(tab, vi_promax, hu_fl) {
  # LEC annotations
  hu_fl$var_id <- paste0("h", seq_len(nrow(hu_fl)))
  tab %>%
    left_join(vi_promax[, 1:5], by = c("var_id" = "UniProt_acc")) %>%
    left_join(hu_fl[, c("UniProt_acc","gene_symbol","var_id")], by = c("antibody" = "var_id")) %>%
    as.data.table()
}

# Key used for replication check
make_qtl_key <- function(dt) {
  dt$direction_of_effect <- ifelse(dt$Beta > 0, "positive", "negative")
  dt$qtl <- paste0(dt$var_id, "_", dt$antibody, "")
  dt$qtl1 <- paste0(dt$var_id, "_", dt$antibody, "_", dt$direction_of_effect)
  dt
}

# -------------------------- LLF block --------------------------
if (!all(c(opt$llf_phe, opt$llf_vi_bin, opt$llf_hu_bin, opt$llf_vi_promax, opt$llf_hu_fl) %in% "") &&
    !any(is.null(c(opt$llf_phe, opt$llf_vi_bin, opt$llf_hu_bin, opt$llf_vi_promax, opt$llf_hu_fl)))) {
  
  message("=== LLF: reading inputs ===")
  llf_phe      <- fread(opt$llf_phe)
  llf_vi_bin   <- fread(opt$llf_vi_bin)     # columns: virus proteins + Subject_Id
  llf_hu_bin   <- fread(opt$llf_hu_bin)     # columns: human FL var_id + Subject_Id
  llf_vi_annot <- fread(opt$llf_vi_promax)
  llf_hu_annot <- fread(opt$llf_hu_fl)
  
  # Merge: phenotype + virus + human
  llf_merged <- llf_phe %>%
    right_join(llf_vi_bin, by = c("Sample_ID"="Subject_Id")) %>%
    right_join(llf_hu_bin, by = c("Sample_ID"="Subject_Id")) %>%
    as.data.table()
  
  # Variables
  hu_cols_llf <- setdiff(colnames(llf_hu_bin), "Subject_Id") #2258
  vi_cols_llf <- setdiff(colnames(llf_vi_bin), "Subject_Id") #2323
  covars_llf  <- c("Age_at_collect", "Gender", "BMI_most_recent", "Race_Group", "Smoking", "Alcohol", "ics_trim_totnum_5y")
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
  
  if ("gene_symbol" %in% names(llf_res)) {
    llf_res$gene_symbol <- ifelse(
      grepl("^ (Unverified)", llf_res$gene_symbol), "", llf_res$gene_symbol)
    llf_res$gene_symbol <- ifelse(
      grepl("^ZNF559-ZNF177,ZNF177", llf_res$gene_symbol), "ZNF177", llf_res$gene_symbol)
  }
  
  # Save all + sig
  out_all_llf <- file.path(opt$out_dir, "llf_hsv_bin_glm_all.tsv")
  out_sig_llf <- file.path(opt$out_dir, "llf_hsv_bin_glm_sig.tsv")
  fwrite(llf_res, out_all_llf, sep = "\t")
  fwrite(llf_res[P.adj < 0.05], out_sig_llf, sep = "\t")
  
  # Report cutoff
  sig_idx <- which(llf_res$P.adj <= 0.05) #0.0001869036
  if (length(sig_idx)) {
    pval_cutoff <- max(llf_res$P[sig_idx], na.rm = TRUE)
    message(sprintf("LLF: FDR<=0.05 p-value cutoff ≈ %.10f", pval_cutoff))
  } else {
    message("LLF: No FDR-significant results at 0.05.")
  }
}

hsv.sig = subset(llf_res, llf_res$P.adj < 0.05) #19608

counts_summarise <- hsv.sig %>% group_by(taxon_species) %>% summarise(total_count = n())
# taxon_species            total_count
# 1 Human alphaherpesvirus 1        2796
# 2 Human alphaherpesvirus 2        2427
# 3 Human alphaherpesvirus 3         254
# 4 Human betaherpesvirus 5        11105
# 5 Human betaherpesvirus 6A         286
# 6 Human betaherpesvirus 6B         324
# 7 Human betaherpesvirus 7           84
# 8 Human gammaherpesvirus 4        2087
# 9 Human gammaherpesvirus 8         245

counts_summarise <- hsv.sig[,c('taxon_species','var_id')] %>% group_by(taxon_species) %>% summarize(unique_count = n_distinct(var_id))
length(unique(hsv.sig$var_id)) #2202
# 1 Human alphaherpesvirus 1          457
# 2 Human alphaherpesvirus 2          321
# 3 Human alphaherpesvirus 3           71
# 4 Human betaherpesvirus 5           668
# 5 Human betaherpesvirus 6A           66
# 6 Human betaherpesvirus 6B           58
# 7 Human betaherpesvirus 7            21
# 8 Human gammaherpesvirus 4          491
# 9 Human gammaherpesvirus 8           49

hsv.sig1= hsv.sig[complete.cases(hsv.sig[ ,c('gene_symbol')]),] #19075
length(unique(hsv.sig1$gene_symbol)) #995
counts_summarise <- hsv.sig1[,c('taxon_species','gene_symbol')] %>% group_by(taxon_species) %>% summarize(unique_count = n_distinct(gene_symbol))
# 1 Human alphaherpesvirus 1          274
# 2 Human alphaherpesvirus 2          298
# 3 Human alphaherpesvirus 3          108
# 4 Human betaherpesvirus 5           467
# 5 Human betaherpesvirus 6A          106
# 6 Human betaherpesvirus 6B           94
# 7 Human betaherpesvirus 7            50
# 8 Human gammaherpesvirus 4          204
# 9 Human gammaherpesvirus 8           76


# -------------------------- LEC block --------------------------
if (!all(c(opt$lec_phe, opt$lec_vi_bin, opt$lec_hu_bin, opt$llf_vi_promax, opt$llf_hu_fl) %in% "") &&
    !any(is.null(c(opt$lec_phe, opt$lec_vi_bin, opt$lec_hu_bin, opt$llf_vi_promax, opt$llf_hu_fl)))) {
  
  message("=== LEC: reading inputs ===")
  lec_phe      <- fread(opt$lec_phe)
  lec_vi_bin   <- fread(opt$lec_vi_bin)     # columns: virus proteins + Sample_Id
  lec_hu_bin   <- fread(opt$lec_hu_bin)     # columns: human FL var_id + Sample_Id

  # Merge: phenotype uses Subject_Id; binaries use Sample_Id
  if (!("Subject_Id" %in% names(lec_phe))) {
    stop("LEC phenotype must contain 'Subject_Id' to join LEC binaries.")
  }
  lec_merged <- lec_phe %>%
    right_join(lec_vi_bin,  by = c("Sample_Id")) %>%
    right_join(lec_hu_bin,  by = c("Sample_Id")) %>%
    as.data.table()

  # Variables
  hu_cols_lec <- setdiff(colnames(lec_hu_bin), "Sample_Id") #1064
  vi_cols_lec <- setdiff(colnames(lec_vi_bin), "Sample_Id") #2342
  covars_lec  <- c("Age_at_collect", "Gender", "BMI_most_recent", "Race_Group", "Smoking", "Alcohol")
  missing_cov <- setdiff(covars_lec, names(lec_merged))
  if (length(missing_cov)) stop("LEC: missing covariates: ", paste(missing_cov, collapse=", "))
  
  message("=== LEC: running GLMs in parallel ===")
  lec_res <- run_glm_parallel(
    df      = lec_merged,
    hu_vars = hu_cols_lec,
    vi_vars = vi_cols_lec,
    covars  = covars_lec,
    n_cores = opt$n_cores
  )
  
  setorder(lec_res, P)
  # Annotate + FDR + direction
  lec_res <- annotate_lec(lec_res, vi_promax = llf_vi_annot, hu_fl = llf_hu_annot)
  lec_res <- make_qtl_key(lec_res)
  
  if ("gene_symbol" %in% names(llf_res)) {
    llf_res$gene_symbol <- ifelse(
      grepl("^ (Unverified)", llf_res$gene_symbol), "", llf_res$gene_symbol)
    llf_res$gene_symbol <- ifelse(
      grepl("^ZNF559-ZNF177,ZNF177", llf_res$gene_symbol), "ZNF177", llf_res$gene_symbol)
  }
  
  # Save all + sig
  out_all_lec <- file.path(opt$out_dir, "lec_hsv_bin_glm_all.tsv")
  out_sig_lec <- file.path(opt$out_dir, "lec_hsv_bin_glm_sig.tsv")
  fwrite(lec_res, out_all_lec, sep = "\t")
  fwrite(lec_res[P < 0.05], out_sig_lec, sep = "\t")

  # ------------------ Replication (LLF -> LEC) ------------------
  if (exists("llf_res")) {
    # LLF significant with gene symbols cleaned like your code
    llf_res = llf_res[complete.cases(llf_res[ ,c('gene_symbol')]),] 
    llf_sig <- copy(llf_res[P.adj < 0.05])

    # LEC significant set
    lec_sig = subset(lec_res, lec_res$qtl%in%hsv.sig1$qtl) #13013 everything included
    lec_sig$P.adj <- p.adjust(lec_sig$P, method = 'fdr')
    
    lec_sig1 = subset(lec_sig, lec_sig$P.adj <0.05) #4784

    # Key already built: qtl = var_id_antibody_direction
    llf_keys <- unique(hsv.sig1$qtl)
    lec_keys <- unique(lec_sig1$qtl)
    common   <- intersect(llf_keys, lec_keys) #4784
    
    rep_llf <- hsv.sig1[qtl %in% common]
    rep_lec <- lec_sig[qtl %in% common]
    
    rep_merge = inner_join(rep_llf, rep_lec, by=c('var_id','antibody','gene_symbol','taxon_genus','taxon_species',
                                                  'UniProt_acc', 'product', 'direction_of_effect', 'qtl'))

    # Save exactly as your filenames
    fwrite(lec_sig, file.path(opt$out_dir, "hsv_bin_fchange_test_lec.tsv"), sep = "\t")
    fwrite(rep_llf, file.path(opt$out_dir, "hsv_bin_fchange_rep_llf.tsv"), sep = "\t")
    fwrite(rep_lec, file.path(opt$out_dir, "hsv_bin_fchange_rep_lec.tsv"), sep = "\t")
    fwrite(rep_merge, file.path(opt$out_dir, "hsv_bin_fchange_rep_both.tsv"), sep = "\t")
    message("Replication complete: N = ", length(common))
  } else {
    message("Replication skipped (LLF results not available in this run).")
  }
}

all(rep_merge$qtl1.x == rep_merge$qtl1.y) #4784 confirm direction
# [1] TRUE


# total associations
counts_total <- lec_sig %>% group_by(taxon_species) %>% summarise(total_count = n())
# 1 Human alphaherpesvirus 1        2000
# 2 Human alphaherpesvirus 2        1221
# 3 Human alphaherpesvirus 3         102
# 4 Human betaherpesvirus 5         8021
# 5 Human betaherpesvirus 6A          83
# 6 Human betaherpesvirus 6B          82
# 7 Human betaherpesvirus 7           19
# 8 Human gammaherpesvirus 4        1426
# 9 Human gammaherpesvirus 8          59

counts_sig <- rep_merge %>% group_by(taxon_species) %>% summarise(total_count = n())
# taxon_species            total_count
# 1 Human alphaherpesvirus 1         648
# 2 Human alphaherpesvirus 2         265
# 3 Human alphaherpesvirus 3           8
# 4 Human betaherpesvirus 5         3575
# 5 Human betaherpesvirus 6A          19
# 6 Human betaherpesvirus 6B          15
# 7 Human gammaherpesvirus 4         254

counts_sum = inner_join(counts_sig, counts_total, by=c('taxon_species'))
counts_sum$pct = counts_sum$total_count.x*100/counts_sum$total_count.y


# total viral peptides
counts_total <- lec_sig[,c('taxon_species','var_id')] %>% group_by(taxon_species) %>% summarize(unique_count = n_distinct(var_id))
# 1 Human alphaherpesvirus 1          417
# 2 Human alphaherpesvirus 2          277
# 3 Human alphaherpesvirus 3           53
# 4 Human betaherpesvirus 5           579
# 5 Human betaherpesvirus 6A           42
# 6 Human betaherpesvirus 6B           42
# 7 Human betaherpesvirus 7             9
# 8 Human gammaherpesvirus 4          434
# 9 Human gammaherpesvirus 8           31

counts_sig <- rep_merge[,c('taxon_species','var_id')] %>% group_by(taxon_species) %>% summarize(unique_count = n_distinct(var_id))
# 1 Human alphaherpesvirus 1          277
# 2 Human alphaherpesvirus 2          171
# 3 Human alphaherpesvirus 3            6
# 4 Human betaherpesvirus 5           383
# 5 Human betaherpesvirus 6A           11
# 6 Human betaherpesvirus 6B            9
# 7 Human gammaherpesvirus 4          164

counts_sum = inner_join(counts_sig, counts_total, by=c('taxon_species'))
counts_sum$pct = counts_sum$unique_count.x*100/counts_sum$unique_count.y

# Gene symbol
counts_total <- lec_sig[,c('taxon_species','gene_symbol')] %>% group_by(taxon_species) %>% summarize(unique_count = n_distinct(gene_symbol))
# 1 Human alphaherpesvirus 1          117
# 2 Human alphaherpesvirus 2          118
# 3 Human alphaherpesvirus 3           29
# 4 Human betaherpesvirus 5           190
# 5 Human betaherpesvirus 6A           35
# 6 Human betaherpesvirus 6B           37
# 7 Human betaherpesvirus 7            13
# 8 Human gammaherpesvirus 4           80
# 9 Human gammaherpesvirus 8           18

counts_sig <- rep_merge[,c('taxon_species','gene_symbol')] %>% group_by(taxon_species) %>% summarize(unique_count = n_distinct(gene_symbol))
# 1 Human alphaherpesvirus 1           17
# 2 Human alphaherpesvirus 2           18
# 3 Human alphaherpesvirus 3            5
# 4 Human betaherpesvirus 5            49
# 5 Human betaherpesvirus 6A            4
# 6 Human betaherpesvirus 6B            4
# 7 Human gammaherpesvirus 4           23

counts_sum = inner_join(counts_sig, counts_total, by=c('taxon_species'))
counts_sum$pct = counts_sum$unique_count.x*100/counts_sum$unique_count.y

length(unique(lec_sig$gene_symbol)) #404
length(unique(lec_sig$var_id)) #1884
length(unique(rep_merge$gene_symbol)) #102
length(unique(rep_merge$var_id)) #1021

message("Done. Results written to: ", normalizePath(opt$out_dir))
