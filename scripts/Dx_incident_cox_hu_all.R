#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(dplyr); library(stringr); library(purrr); library(tidyr);
  library(ggplot2); library(ggrepel); library(broom); library(survival); library(coxphf);library(meta)   # Firth Cox for rare events
  source("R/utils_common.R"); source("R/utils_qtl.R")
})

# ----------------------------- CLI -----------------------------
option_list <- list(
  make_option("--llf_pheno",      type="character", help="llf_1290_phe.tsv"),
  make_option("--llf_vi_promax",   type="character", help="IB1007_VirSIGHT_Promax_Hits_Fold-Over-Background.csv"),
  make_option("--llf_hu_fl",   type="character", help="IB1007_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv"),
  make_option("--llf_vi_bin",   type="character", help="hsv_promax_bin_MGBB-LLF.tsv"),
  make_option("--llf_hu_bin",   type="character", help="human_fl_bin_MGBB-LLF.tsv"),
  make_option("--rep_llf",   type="character", help="hsv_bin_fchange_rep_llf.tsv"),
  make_option("--llf_dx_count",   type="character", help="llf_case_counts_filtered.csv"),
  make_option("--lec_pheno",   type="character", help="lec_763_phe.tsv"),
  make_option("--lec_vi_bin",   type="character", help="hsv_promax_bin_MGBB-LEC.tsv"),
  make_option("--lec_hu_bin",   type="character", help="human_fl_bin_MGBB-LEC.tsv"),
  make_option("--abc_vi_promax",   type="character", help="IB1021_VirSIGHT_Promax_Hits_Fold-Over-Background.csv"),
  make_option("--abc_hu_fl",   type="character", help="IB1021_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv"),
  make_option("--leo_vi_promax",   type="character", help="IB1189_VirSIGHT_Promax_Hits_Fold-Over-Background_Old-Format-HSV.tsv"),
  make_option("--leo_hu_fl",   type="character", help="IB1189_HuSIGHT_FullLength_Hits_Fold-Over-Background.tsv"),
  make_option("--rep_lec",   type="character", help="hsv_bin_fchange_rep_lec.tsv"),
  make_option("--lec_dx_count",   type="character", help="lec_case_counts_filtered.csv"),
  make_option("--sig_ab",   type="character", help="lec_test_auc_summary_sig.csv"),
  make_option("--out_dir",  type="character", default="results/Disease", help="Output dir")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

# Update function to accept cohort name and data explicitly
run_cox_human <- function(ii, cohort_name, input_dt) {
  
  dx <- dx_list[ii]
  time_col  <- c(paste0("time_to_", gsub('_flag', '', dx)))
  
  print(paste("Processing:", dx,"", time_col,"| Cohort:", cohort_name))
  
  ## ---------- 1. Basic Counts ----------------------------------------------
  # Check if disease and time columns exist
  if (!dx %in% names(input_dt)) {
    warning(paste("Disease", dx, "not found in cohort", cohort_name))
    return(NULL)
  }
  if (!time_col %in% names(input_dt)) {
    stop(paste("Time column", time_col, "not found in data! Please update 'time_col' variable."))
  }
  
  case_idx    <- input_dt[[dx]] == 1
  control_idx <- input_dt[[dx]] == 0
  n_case      <- sum(case_idx,    na.rm = TRUE)
  n_control   <- sum(control_idx, na.rm = TRUE)
  
  # Calculate positive counts (Protein Presence)
  case_pos_vec    <- colSums(input_dt[case_idx,    ..hu_var] > 1, na.rm = TRUE)
  control_pos_vec <- colSums(input_dt[control_idx, ..hu_var] > 1, na.rm = TRUE)
  
  if (cohort_name=='llf'){
    covars <- "Age_at_collect + Gender + BMI_most_recent + Race_Group + Smoking + Alcohol + ics_trim_totnum_5y + cci"
  } else {covars <- "Age_at_collect + Gender + BMI_most_recent + Race_Group + Smoking + Alcohol + cci"}
  
  ## ---------- 2. Per-Protein Loop -----------------------------------------
  cox_results <- lapply(seq_along(hu_var), function(j) {
    
    var <- hu_var[j]
    
    # --- Check for Sparsity (< 5 in any cell) ---
    # Even in Cox, we check 2x2 sparsity as a proxy for stability
    c_pos <- case_pos_vec[j]
    c_neg <- n_case - c_pos
    u_pos <- control_pos_vec[j]
    u_neg <- n_control - u_pos
    
    is_sparse <- min(c_pos, c_neg, u_pos, u_neg) < 5
    
    # --- Create Formula ---
    # Surv(time, event) ~ var + covars
    fml_str <- paste0("Surv(", time_col, ", ", dx, ") ~ ", var, " + ", covars)
    fml     <- as.formula(fml_str)
    
    # Initialize fields
    beta <- se <- z <- p <- AIC_val <- NA_real_
    converged <- FALSE
    model_method <- "coxph" 
    
    # --- Fit Model ---
    res <- tryCatch({
      if (is_sparse) {
        model_method <- "coxphf"
        # coxphf: Firth's Penalized Likelihood for Cox
        # pl=FALSE for Wald (faster), pl=TRUE for Profile Likelihood (slower/accurate)
        
        covar_list <- unlist(strsplit(covars, " \\+ "))
        input_dt.1 = na.omit(subset(input_dt, select=c(dx, time_col, var, covar_list)))
        
        coxphf(fml, data = input_dt.1, pl = FALSE)
      } else {
        # Standard Cox
        coxph(fml, data = input_dt)
      }
    }, error = function(e) return(NULL))
    
    
    # --- Extract Results ---
    if (!is.null(res)) {
      
      if (is_sparse) {
        # --- Extraction for Firth Cox (coxphf) ---
        converged <- TRUE 
        
        # coxphf stores coefficients in a simple vector
        coef_idx <- which(names(res$coefficients) == var)
        
        if (length(coef_idx) > 0) {
          beta <- res$coefficients[coef_idx]
          # Extract SE from the variance-covariance matrix
          se   <- sqrt(diag(res$var))[coef_idx]
          p    <- res$prob[coef_idx]
          z    <- beta / se 
        }
        # AIC is calculated differently for penalized models, often omitted or approximated
        AIC_val <- -2 * res$loglik[2] + 2 * res$df
        
      } else {
        # --- Extraction for Standard Cox (coxph) ---
        # Checks for convergence
        if (!is.null(res$iter)) converged <- TRUE 
        
        # Summary structure for coxph:
        # Columns: 1=coef, 2=exp(coef), 3=se(coef), 4=z, 5=p
        cf <- summary(res)$coefficients
        
        if (var %in% rownames(cf)) {
          beta <- cf[var, 1] 
          se   <- cf[var, 3] 
          z    <- cf[var, 4] 
          p    <- cf[var, 5]
        } else if (nrow(cf) >= 1) {
          # Fallback if name exact match fails but var is first predictor
          beta <- cf[1,1]; se <- cf[1,3]; z <- cf[1,4]; p <- cf[1,5]
        }
        
        AIC_val <- extractAIC(res)[2]
      }
    }
    
    data.table(
      Beta        = beta,
      HR          = exp(beta), # Added Hazard Ratio
      SE          = se,
      Z           = z,
      P           = p,
      n_case      = n_case,
      n_control   = n_control,
      case_pos    = c_pos,
      case_neg    = c_neg,
      control_pos = u_pos,
      control_neg = u_neg,
      AIC         = AIC_val,
      converged   = converged,
      var_id      = var,
      is_sparse   = is_sparse,
      method      = model_method,
      cohort      = cohort_name
    )
  })
  
  out <- rbindlist(cox_results, use.names = TRUE, fill = TRUE)
  out[, disease := dx]
  setorder(out, P)
  
  # Filename updated to 'cox'
  outfile <- file.path(opt$out_dir, paste0("bin_fchange_hu_", dx, "_cox_", cohort_name,".csv"))
  fwrite(out, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE)
}


# read all result files from a pattern and annotate disease
read_results <- function(pattern) {
  fs <- list.files(path = opt$out_dir, pattern = pattern, full.names = TRUE)
  if (!length(fs)) return(data.table())
  rbindlist(lapply(fs, function(f) {
    dt <- fread(f)
    # infer disease/sex from filename if missing
    if (!"disease" %in% names(dt)) {
      m <- regmatches(f, regexec("_hu_(.+?)_", f))[[1]]
      if (length(m) >= 2) dt[, disease := m[2]]
    }
    if (!"cohort" %in% names(dt)) {
      if (grepl("_llf\\.csv", f)) dt[, cohort := "llf"]
      if (grepl("_lec\\.csv", f)) dt[, cohort := "lec"]
      if (grepl("_all\\.csv", f)) dt[, cohort := "all"]
    }
    dt[, file := basename(f)]
    dt
  }), use.names = TRUE, fill = TRUE)
}

# ----------------------------- IO in LLF ------------------------------
# Phenotype table and covariate processing in LLF
llf_phe <- fread(opt$llf_pheno)
llf_phe[, Subject_Id := as.character(Subject_Id)]
llf_phe1 <- build_disease_combine(llf_phe)

# AUC
sig_ab <- fread(opt$sig_ab)

# MIPSA files
llf_vi_promax <- fread(opt$llf_vi_promax)
llf_hu_fl <- fread(opt$llf_hu_fl)
if (!"var_id" %in% names(llf_hu_fl)) llf_hu_fl[, var_id := paste0("h", seq_len(.N))]

# Select subjects and Merge
sample_lists = llf_phe1$Sample_ID
llf.fl = subset(llf_hu_fl, select = sample_lists) #1290 subjects
llf.fl[llf.fl==1] <- 0

llf.fl1 = as.data.frame(t(llf.fl))
colnames(llf.fl1) = llf_hu_fl$var_id
llf.fl1$Sample_ID = rownames(llf.fl1)

merge_dt_raw.llf <- llf_phe1 %>%
  right_join(llf.fl1, by = "Sample_ID") %>%
  as.data.table()

# ----------------------------- Run COX -----------------------------
# Human FL predictors (all antibodies)
rep_llf = fread(opt$rep_llf)
hu_var <- unique(rep_llf$antibody) #117

# ----------------------------- Disease cols -----------------------------
dx_cols <- grep("_combine$", colnames(merge_dt_raw.llf), value = TRUE)
if (length(dx_cols) == 0) stop("No *_combine disease columns found. Did build_disease_combine() run correctly?")

# Keep diseases with at least 13 positives
dx_sums <- colSums(merge_dt_raw.llf[, ..dx_cols], na.rm = TRUE)
dx_list <- names(dx_sums[dx_sums >= 12.9])
ex_list <- c(default="AIDS","Depression","HepatitisB","HepatitisC","Migraine","Opioid_use_disorder","Headache","Asthma","Chronic_viral_hepatitis")
dx_list <- setdiff(dx_list, paste0(ex_list,"_combine")) # 62 diseases
dx_list <- gsub("_combine", "", dx_list)
dx_list <- paste0(dx_list,"_flag")

lapply(1:length(dx_list), run_cox_human, cohort_name = "llf", input_dt = merge_dt_raw.llf)

# ----------------------------- IO in LEC------------------------------
# Phenotype table and covariate processing in LEC
lec_phe <- fread(opt$lec_pheno)
lec_phe1 <- build_disease_combine(lec_phe)

# MIPSA files
abc_vi_promax <- fread(opt$abc_vi_promax)
abc_hu_fl <- fread(opt$abc_hu_fl)
if (!"var_id" %in% names(abc_hu_fl)) abc_hu_fl[, var_id := paste0("h", .I)]

leo_vi_promax <- fread(opt$leo_vi_promax)
leo_hu_fl <- fread(opt$leo_hu_fl)

abc_hu_fl1 = subset(abc_hu_fl, abc_hu_fl$pep_id%in%leo_hu_fl$pep_id)
all(leo_hu_fl$pep_id==abc_hu_fl1$pep_id)
# [1] TRUE

lec_hu_fl = cbind(abc_hu_fl1, leo_hu_fl)

sample_lists = lec_phe1$Sample_Id
lec.fl <- subset(lec_hu_fl, select = sample_lists)
lec.fl[lec.fl==1] <- 0

lec.fl1 = as.data.frame(t(lec.fl))
colnames(lec.fl1) = lec_hu_fl$var_id
lec.fl1$Sample_Id = rownames(lec.fl1)

merge_dt_raw.lec <- lec_phe1 %>%
  right_join(lec.fl1, by = "Sample_Id") %>%
  as.data.table()

# ----------------------------- Disease cols -----------------------------
# Select the same diseases in LLF tested 
dx_test <- colSums(merge_dt_raw.lec[, ..dx_list], na.rm = TRUE)
length(dx_test)==length(dx_list)

lapply(1:length(dx_list), run_cox_human, cohort_name = "lec", input_dt = merge_dt_raw.lec)

# ----------------------------- IO in Pooled LLF-LEC-----------------------------
common_names = intersect(colnames(merge_dt_raw.lec), colnames(merge_dt_raw.llf))

merge_dt_raw.all = rbind(merge_dt_raw.llf[, ..common_names], merge_dt_raw.lec[, ..common_names]) #2053

dx_test <- colSums(merge_dt_raw.all[, ..dx_list], na.rm = TRUE)
length(dx_test)==length(dx_list)

lapply(1:length(dx_list), run_cox_human, cohort_name = "all", input_dt = merge_dt_raw.all)

# ----------------------------- Collect Data: Every cohorts -----------------------------

hu.dx_all <- read_results("^bin_fchange_hu_.*_(llf|lec|all)\\.csv$")
hu.dx_all <- hu.dx_all[,-c('file')]
hu.dx_all$disease <- gsub('_flag', '', hu.dx_all$disease) 
hu.dx_all = inner_join(llf_hu_fl[,c('var_id', 'gene_symbol','UniProt_acc', 'product')], hu.dx_all, by=c('var_id'))
hu.dx_all$gene_symbol <- gsub(" \\s*\\([^\\)]+\\)", "",  hu.dx_all$gene_symbol)
hu.dx_all$gene_symbol <- gsub("ZNF559-ZNF177,ZNF177", "ZNF177",  hu.dx_all$gene_symbol)

## Reshape data to wide format 
res_wide <- dcast(hu.dx_all, disease + var_id + gene_symbol + UniProt_acc + product ~ cohort, 
                  value.var = c("HR","Beta","SE","Z", "P","n_case",'n_control',"case_pos", "control_pos"))

dx_lookup <- stack(dx_categories)
names(dx_lookup) <- c("Disease", "Category")
res_wide = left_join(res_wide, dx_lookup, by=c('disease'='Disease'))

res_wide$disease = gsub('_', ' ', res_wide$disease)
res_wide$Category = gsub('_', ' ', res_wide$Category)

res_wide$p.adj_lec = p.adjust(res_wide$P_lec, method='fdr')
res_wide$p.adj_all = p.adjust(res_wide$P_all, method='fdr')
res_wide$p.adj_llf = p.adjust(res_wide$P_llf, method='fdr')

fwrite(res_wide,
       file.path(opt$out_dir, "rep_hu_dx_all_inc_coxphf_ann.tsv"),
       sep = "\t", quote = FALSE)
