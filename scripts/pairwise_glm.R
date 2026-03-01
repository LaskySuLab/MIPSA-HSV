#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(parallel)
  source("R/utils_common.R"); source("R/utils_qtl.R")
})

# ---------------------------- CLI ----------------------------
opt_list <- list(
  # LLF inputs
  make_option("--llf_phe", type="character",
              help="llf_1289_phe.tsv"),
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
    right_join(llf_vi_bin, by = c("Sample_Id"="Subject_Id")) %>%
    right_join(llf_hu_bin, by = c("Sample_Id"="Subject_Id")) %>%
    as.data.table()
  
  # Variables
  hu_cols_llf <- setdiff(colnames(llf_hu_bin), "Subject_Id") #2257
  vi_cols_llf <- setdiff(colnames(llf_vi_bin), "Subject_Id") #2322
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
      grepl("\\s*\\(Unverified\\)", llf_res$gene_symbol), "", llf_res$gene_symbol)
  }
  
  # Save all + sig
  out_all_llf <- file.path(opt$out_dir, "llf_hsv_bin_glm_all.tsv")
  out_sig_llf <- file.path(opt$out_dir, "llf_hsv_bin_glm_sig.tsv")
  fwrite(llf_res, out_all_llf, na='', sep = "\t")
  fwrite(llf_res[P.adj < 0.05], out_sig_llf, na='', sep = "\t")
  
  # Report cutoff
  sig_idx <- which(llf_res$P.adj <= 0.05)
  if (length(sig_idx)) {
    pval_cutoff <- max(llf_res$P[sig_idx], na.rm = TRUE)
    message(sprintf("LLF: FDR<=0.05 p-value cutoff ≈ %.10f", pval_cutoff))
  } else {
    message("LLF: No FDR-significant results at 0.05.")
  }
}

hsv.sig = subset(llf_res, llf_res$P.adj < 0.05) #19812

counts_summarise <- hsv.sig %>% group_by(taxon_species) %>% summarise(total_count = n())
# taxon_species            total_count
# 1 Human alphaherpesvirus 1        2845
# 2 Human alphaherpesvirus 2        2369
# 3 Human alphaherpesvirus 3         258
# 4 Human betaherpesvirus 5        11307
# 5 Human betaherpesvirus 6A         293
# 6 Human betaherpesvirus 6B         331
# 7 Human betaherpesvirus 7           83
# 8 Human gammaherpesvirus 4        2062
# 9 Human gammaherpesvirus 8         264

counts_summarise <- hsv.sig[,c('taxon_species','var_id')] %>% group_by(taxon_species) %>% summarize(unique_count = n_distinct(var_id))
length(unique(hsv.sig$var_id)) #2194
# 1 Human alphaherpesvirus 1          458
# 2 Human alphaherpesvirus 2          322
# 3 Human alphaherpesvirus 3           76
# 4 Human betaherpesvirus 5           667
# 5 Human betaherpesvirus 6A           65
# 6 Human betaherpesvirus 6B           58
# 7 Human betaherpesvirus 7            20
# 8 Human gammaherpesvirus 4          479
# 9 Human gammaherpesvirus 8           49

hsv.sig1= hsv.sig[complete.cases(hsv.sig[ ,c('gene_symbol')]),] #19270
length(unique(hsv.sig1$gene_symbol)) #990
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
    right_join(lec_vi_bin,  by = c("Sample_Id"="Subject_Id")) %>%
    right_join(lec_hu_bin,  by = c("Sample_Id"="Subject_Id")) %>%
    as.data.table()

  # Variables
  hu_cols_lec <- setdiff(colnames(lec_hu_bin), "Sample_Id") #1065
  vi_cols_lec <- setdiff(colnames(lec_vi_bin), "Sample_Id") #2322
  covars_lec  <- c("Age_at_collect", "Gender", "BMI_most_recent", "Race_Group", "Smoking", "Alcohol", "cohort")
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
      grepl("\\s*\\(Unverified\\)", llf_res$gene_symbol), "", llf_res$gene_symbol)
  }
  
  # Save all + sig
  out_all_lec <- file.path(opt$out_dir, "lec_hsv_bin_glm_all.tsv")
  out_sig_lec <- file.path(opt$out_dir, "lec_hsv_bin_glm_sig.tsv")
  fwrite(lec_res, out_all_lec, na='', sep = "\t")
  fwrite(lec_res[P < 0.05], out_sig_lec, na='', sep = "\t")

  # ------------------ Replication (LLF -> LEC) ------------------
  if (exists("llf_res")) {
    # LLF significant with gene symbols cleaned like your code
    llf_res1 = llf_res[complete.cases(llf_res[ ,c('gene_symbol')]),] 
    llf_sig <- copy(llf_res1[P.adj < 0.05])

    # LEC significant set
    lec_sig = subset(lec_res, lec_res$qtl%in%hsv.sig1$qtl) #13132 everything included
    lec_sig$P.adj <- p.adjust(lec_sig$P, method = 'fdr')
    
    lec_sig1 = subset(lec_sig, lec_sig$P.adj <0.05) #3948

    # Key already built: qtl = var_id_antibody_direction
    llf_keys <- unique(hsv.sig1$qtl1)
    lec_keys <- unique(lec_sig1$qtl1)
    common   <- intersect(llf_keys, lec_keys) #3943
    
    rep_llf <- hsv.sig1[qtl1 %in% common]
    rep_lec <- lec_sig[qtl1 %in% common]
    
    rep_merge = inner_join(rep_llf, rep_lec, by=c('var_id','antibody','gene_symbol','taxon_genus','taxon_species',
                                                  'UniProt_acc', 'product', 'direction_of_effect', 'qtl','qtl1'))

    # Save exactly as your filenames
    fwrite(lec_sig, file.path(opt$out_dir, "hsv_bin_fchange_test_lec.tsv"), na='', sep = "\t")
    fwrite(rep_llf, file.path(opt$out_dir, "hsv_bin_fchange_rep_llf.tsv"), na='', sep = "\t")
    fwrite(rep_lec, file.path(opt$out_dir, "hsv_bin_fchange_rep_lec.tsv"), na='', sep = "\t")
    fwrite(rep_merge, file.path(opt$out_dir, "hsv_bin_fchange_rep_both.tsv"), na='', sep = "\t")
    message("Replication complete: N = ", length(common))
  } else {
    message("Replication skipped (LLF results not available in this run).")
  }
}

# total associations
counts_total <- lec_sig %>% group_by(taxon_species) %>% summarise(total_count = n())
# 1 Human alphaherpesvirus 1        2049
# 2 Human alphaherpesvirus 2        1180
# 3 Human alphaherpesvirus 3         109
# 4 Human betaherpesvirus 5         8152
# 5 Human betaherpesvirus 6A          86
# 6 Human betaherpesvirus 6B          85
# 7 Human betaherpesvirus 7           18
# 8 Human gammaherpesvirus 4        1386
# 9 Human gammaherpesvirus 8          67

counts_sig <- rep_merge %>% group_by(taxon_species) %>% summarise(total_count = n())
counts_sum = inner_join(counts_sig, counts_total, by=c('taxon_species'))
counts_sum$pct = counts_sum$total_count.x*100/counts_sum$total_count.y
# 1 Human alphaherpesvirus 1           538          2049 26.3 
# 2 Human alphaherpesvirus 2           180          1180 15.3 
# 3 Human alphaherpesvirus 3             5           109  4.59
# 4 Human betaherpesvirus 5           3050          8152 37.4 
# 5 Human betaherpesvirus 6A             9            86 10.5 
# 6 Human betaherpesvirus 6B             7            85  8.24
# 7 Human gammaherpesvirus 4           153          1386 11.0 
sum(counts_sum$total_count.x)
# [1] 3942
sum(counts_sum$total_count.y) 
# [1] 13047

# total viral peptides
counts_total <- lec_sig[,c('taxon_species','var_id')] %>% group_by(taxon_species) %>% summarize(unique_count = n_distinct(var_id))
# 1 Human alphaherpesvirus 1          421
# 2 Human alphaherpesvirus 2          276
# 3 Human alphaherpesvirus 3           61
# 4 Human betaherpesvirus 5           585
# 5 Human betaherpesvirus 6A           41
# 6 Human betaherpesvirus 6B           42
# 7 Human betaherpesvirus 7             9
# 8 Human gammaherpesvirus 4          417
# 9 Human gammaherpesvirus 8           37

counts_sig <- rep_merge[,c('taxon_species','var_id')] %>% group_by(taxon_species) %>% summarize(unique_count = n_distinct(var_id))
counts_sum = inner_join(counts_sig, counts_total, by=c('taxon_species'))
counts_sum$pct = counts_sum$unique_count.x*100/counts_sum$unique_count.y
# 1 Human alphaherpesvirus 1            258            421 61.3 
# 2 Human alphaherpesvirus 2            133            276 48.2 
# 3 Human alphaherpesvirus 3              5             61  8.20
# 4 Human betaherpesvirus 5             377            585 64.4 
# 5 Human betaherpesvirus 6A              6             41 14.6 
# 6 Human betaherpesvirus 6B              4             42  9.52
# 7 Human gammaherpesvirus 4             98            417 23.5 
sum(counts_sum$unique_count.x)
# [1] 881
sum(counts_sum$unique_count.y) 
# [1] 1843

# Gene symbol
counts_total <- lec_sig[,c('taxon_species','gene_symbol')] %>% group_by(taxon_species) %>% summarize(unique_count = n_distinct(gene_symbol))
# 1 Human alphaherpesvirus 1          123
# 2 Human alphaherpesvirus 2          116
# 3 Human alphaherpesvirus 3           30
# 4 Human betaherpesvirus 5           183
# 5 Human betaherpesvirus 6A           36
# 6 Human betaherpesvirus 6B           38
# 7 Human betaherpesvirus 7            13
# 8 Human gammaherpesvirus 4           74
# 9 Human gammaherpesvirus 8           19

counts_sig <- rep_merge[,c('taxon_species','gene_symbol')] %>% group_by(taxon_species) %>% summarize(unique_count = n_distinct(gene_symbol))
counts_sum = inner_join(counts_sig, counts_total, by=c('taxon_species'))
counts_sum$pct = counts_sum$unique_count.x*100/counts_sum$unique_count.y
# 1 Human alphaherpesvirus 1             14            123 11.4 
# 2 Human alphaherpesvirus 2             12            116 10.3 
# 3 Human alphaherpesvirus 3              2             30  6.67
# 4 Human betaherpesvirus 5              43            183 23.5 
# 5 Human betaherpesvirus 6A              2             36  5.56
# 6 Human betaherpesvirus 6B              1             38  2.63
# 7 Human gammaherpesvirus 4             19             74 25.7 

sum(counts_sum$unique_count.x)
# [1] 93
sum(counts_sum$unique_count.y) 
# [1] 600

length(unique(lec_sig$gene_symbol)) #400
length(unique(lec_sig$var_id)) #1889
length(unique(rep_merge$gene_symbol)) #80
length(unique(rep_merge$var_id)) #881

message("Done. Results written to: ", normalizePath(opt$out_dir))
