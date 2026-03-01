#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse); library(data.table); library(dplyr); library(stringr); library(purrr); library(tidyr);
  library(ggplot2); library(ggrepel); library(broom); library(meta); library(logistf) # Required for Firth regression
  source("R/utils_common.R"); source("R/utils_qtl.R")
})

# ----------------------------- CLI -----------------------------
option_list <- list(
  make_option("--phe",      type="character", help="llf_1289_phe.tsv"),
  make_option("--vi_bin",   type="character", help="hsv_proteins_binary.txt"),
  make_option("--hu_bin",   type="character", help="human_fl_binary.txt"),
  make_option("--dx_count", type="character",
              help="a directory containing prevalence/incidence files"),
  make_option("--out_dir",  type="character", default="results/Dx_pre", help="Output dir")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

# Function for parallel processing for human antibodies
run_glm_virus <- function(ii, cohort_name, input_dt) {
  
  dx <- dx_list[ii]
  print(paste("Processing:", dx, "| Cohort:", cohort_name))
  
  ## ---------- 1. Basic Counts (Using input_dt) ----------------------------
  # Check if disease exists in this specific dataset
  if (!dx %in% names(input_dt)) {
    warning(paste("Disease", dx, "not found in cohort", cohort_name))
    return(NULL)
  }
  
  case_idx    <- input_dt[[dx]] == 1
  control_idx <- input_dt[[dx]] == 0
  n_case      <- sum(case_idx,    na.rm = TRUE)
  n_control   <- sum(control_idx, na.rm = TRUE)
  
  # Calculate positive counts
  case_pos_vec    <- colSums(input_dt[case_idx,    ..vi_var] > 1, na.rm = TRUE)
  control_pos_vec <- colSums(input_dt[control_idx, ..vi_var] > 1, na.rm = TRUE)
  
  if (cohort_name=='llf'){
    covars <- "Age_at_collect + Gender + BMI_most_recent + Race_Group + Smoking + Alcohol + ics_trim_totnum_5y + cci"
  } else {covars <- "Age_at_collect + Gender + BMI_most_recent + Race_Group + Smoking + Alcohol + cci + cohort"}
  
  ## ---------- 2. Per-Protein Loop -----------------------------------------
  glm_results <- lapply(seq_along(vi_var), function(j) {
    
    var <- vi_var[j]
    
    # --- Check for Sparsity (< 5 in any cell) ---
    c_pos <- case_pos_vec[j]
    c_neg <- n_case - c_pos
    u_pos <- control_pos_vec[j]
    u_neg <- n_control - u_pos
    
    is_sparse <- min(c_pos, c_neg, u_pos, u_neg) < 5
    
    fml <- as.formula(paste0(dx, " ~ ", var, " + ", covars))
    
    # Initialize fields
    beta <- se <- z <- p <- AIC_val <- maxCook <- NA_real_
    converged <- FALSE
    model_method <- "glm" 
    
    # --- Fit Model ---
    res <- tryCatch({
      if (is_sparse) {
        model_method <- "firth"
        logistf(fml, data = input_dt, pl = FALSE) 
      } else {
        glm(fml, family = "binomial", data = input_dt)
      }
    }, error = function(e) return(NULL))
    
    
    # --- Extract Results ---
    if (!is.null(res)) {
      if (is_sparse) {
        # Extraction for Firth
        converged <- TRUE 
        coef_idx <- which(names(res$coefficients) == var)
        if (length(coef_idx) > 0) {
          beta <- res$coefficients[coef_idx]
          se   <- sqrt(diag(vcov(res)))[coef_idx]
          p    <- res$prob[coef_idx]
          z    <- beta / se 
        }
        AIC_val <- -2 * res$loglik[2] + 2 * res$df
        maxCook <- NA_real_
        
      } else {
        # Extraction for Standard GLM
        converged <- isTRUE(res$converged)
        cf <- summary(res)$coefficients
        
        # Robust row matching
        if (var %in% rownames(cf)) {
          beta <- cf[var, 1]; se <- cf[var, 2]; z <- cf[var, 3]; p <- cf[var, 4]
        } else if (nrow(cf) >= 2) {
          beta <- cf[2,1]; se <- cf[2,2]; z <- cf[2,3]; p <- cf[2,4]
        }
        
        AIC_val <- broom::glance(res)$AIC
        maxCook <- suppressWarnings(max(cooks.distance(res)))
      }
    }
    
    data.table(
      Beta        = beta,
      OR          = exp(beta), # Added Odd Ratio
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
      maxCook     = maxCook,
      converged   = converged,
      var_id      = var,
      is_sparse   = is_sparse,
      method      = model_method,
      cohort      = cohort_name
    )
  })
  
  out <- rbindlist(glm_results, use.names = TRUE, fill = TRUE)
  out[, disease := dx]
  setorder(out, P)
  
  # Filename includes cohort_name
  outfile <- file.path(opt$out_dir, paste0("bin_fchange_vi_", dx, "_glm_", cohort_name, ".csv"))
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
      m <- regmatches(f, regexec("_vi_(.+?)_", f))[[1]]
      if (length(m) >= 2) dt[, disease := m[2]]
    }
    if (!"cohort" %in% names(dt)) {
      if (grepl("_llf\\.csv", f))   dt[, cohort := "llf"]
      if (grepl("_ale\\.csv", f)) dt[, cohort := "lec"]
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
sample_lists = llf_phe1$Sample_Id
llf.vi = subset(llf_vi_promax, select = sample_lists) #1289 subjects
llf.vi[llf.vi==1] <- 0

llf.vi1 = as.data.frame(t(llf.vi))
colnames(llf.vi1) = llf_vi_promax$UniProt_acc
llf.vi1$Sample_Id = rownames(llf.vi1)

merge_dt_raw.llf <- llf_phe1 %>%
  right_join(llf.vi1, by = "Sample_Id") %>%
  as.data.table()

# ----------------------------- Run GLMs -----------------------------
# Human FL predictors (all antibodies)
rep_llf = fread(opt$rep_llf)
vi_var <- unique(rep_llf$var_id) #881

# ----------------------------- Disease cols -----------------------------
dx_cols <- grep("_combine$", colnames(merge_dt_raw.llf), value = TRUE)
if (length(dx_cols) == 0) stop("No *_combine disease columns found. Did build_disease_combine() run correctly?")

# Keep diseases with at least 13 positives
dx_sums <- colSums(merge_dt_raw.llf[, ..dx_cols], na.rm = TRUE)
dx_list <- names(dx_sums[dx_sums >= 12.89])
ex_list <- c(default="Depression","Opioid_use_disorder","Alcoholism","toothache","sciatica","Psychosis","loose_tooth",
             "Asthma","Headache","Migraine","Post_Traumatic_Stress_Disorder","Bipolar_Disorder","Female_Infertility","Male_Infertility",
             "Chronic viral hepatitis","Schizophrenia","Chronic_viral_hepatitis","Other_chronic_heptitis","Alcohol_liver_disease","Hiatus_Hernia",
             "Diverticular","Cataract","Presbyopia","Mouth_ulcer")
dx_list <- setdiff(dx_list, paste0(ex_list,"_combine")) # 73 diseases
dx_list <- gsub("_combine", "", dx_list)
dx_list <- paste0(dx_list,"_at_collect")

lapply(1:length(dx_list), run_glm_virus, cohort_name = "llf", input_dt = merge_dt_raw.llf)

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

all(leo_vi_promax$UniProt_acc==abc_vi_promax$UniProt_acc)
# [1] TRUE
# Raw data
lec_vi_promax = cbind(abc_vi_promax, leo_vi_promax)

sample_lists = lec_phe1$Sample_Id
lec.vi <- subset(lec_vi_promax, select = sample_lists)
lec.vi[lec.vi==1] <- 0

lec.vi1 = as.data.frame(t(lec.vi))
colnames(lec.vi1) = lec_vi_promax$UniProt_acc
lec.vi1$Sample_Id = rownames(lec.vi1)

merge_dt_raw.lec <- lec_phe1 %>%
  inner_join(lec.vi1, by = "Sample_Id") %>%
  as.data.table()

# ----------------------------- Disease cols -----------------------------
# Select the same diseases in LLF tested 
dx_test <- colSums(merge_dt_raw.lec[, ..dx_list], na.rm = TRUE)
length(dx_test)==length(dx_list)

lapply(1:length(dx_list), run_glm_virus, cohort_name = "lec", input_dt = merge_dt_raw.lec)

# ----------------------------- IO in Pooled LLF-LEC------------------------------
common_names = intersect(colnames(merge_dt_raw.lec), colnames(merge_dt_raw.llf))

merge_dt_raw.all = rbind(merge_dt_raw.llf[, ..common_names], merge_dt_raw.lec[, ..common_names]) #2053

dx_test <- colSums(merge_dt_raw.all[, ..dx_list], na.rm = TRUE)
length(dx_test)==length(dx_list)

lapply(1:length(dx_list), run_glm_virus, cohort_name = "all", input_dt = merge_dt_raw.all)

# ----------------------------- Collect Data: Every cohorts -----------------------------

vi.dx_all <- read_results("^bin_fchange_vi_.*_(llf|lec|all)\\.csv$")
vi.dx_all <- vi.dx_all[,-c('file')]
vi.dx_all$disease <- gsub('_at_collect', '', vi.dx_all$disease) 
vi.dx_all = inner_join(llf_vi_promax[,c('UniProt_acc','taxon_species','taxon_genus', 'product')], vi.dx_all, by=c('UniProt_acc'='var_id'))

## Reshape data to wide format 
res_wide <- dcast(vi.dx_all, disease + UniProt_acc + taxon_genus + taxon_species + product ~ cohort, 
                  value.var = c("OR","Beta","SE","Z", "P","n_case",'n_control',"case_pos", "control_pos"))

dx_lookup <- stack(dx_categories)
names(dx_lookup) <- c("Disease", "Category")
res_wide = left_join(res_wide, dx_lookup, by=c('disease'='Disease'))

res_wide$disease = gsub('_', ' ', res_wide$disease)
res_wide$Category = gsub('_', ' ', res_wide$Category)

res_wide$p.adj_lec = p.adjust(res_wide$P_lec, method='fdr')
res_wide$p.adj_all = p.adjust(res_wide$P_all, method='fdr')
res_wide$p.adj_llf = p.adjust(res_wide$P_llf, method='fdr')

fwrite(res_wide,
       file.path(opt$out_dir, "rep_vi_dx_all_pre_glmf_ann.tsv"),
       sep = "\t", quote = FALSE)
fwrite(subset(res_wide, res_wide$P_all <0.05) ,
       file.path(opt$out_dir, "rep_vi_dx_all_pre_glmf_ann_sig.tsv"),
       sep = "\t", quote = FALSE)

res_wide = fread(file.path(opt$out_dir, "rep_vi_dx_all_pre_glmf_ann.tsv"), head=T)

# Classify Significance
res_wide[, Sig_Status := fcase(
  P_lec < 0.05 & P_llf < 0.05, "Both Significant",
  P_lec < 0.05,                "MGBB-LEC Only",
  P_llf < 0.05,                "MGBB-LLF Only",
  default =                    "Not Significant"
)]

res_wide$taxon_species <- species_label_map(res_wide$taxon_species)

# Plot
cutoff_p <- {
  sig_idx <- which(res_wide$p.adj_all < 0.05)
  if (length(sig_idx)) max(res_wide$P_all[sig_idx], na.rm = TRUE) else NULL
} #2.821267e-06

p.all <- ggplot(res_wide, aes(x = disease, y = -log10(P_all))) +
  facet_grid(. ~ Category, scales = "free_x") +
  geom_point(aes(color = OR_all), size = 2) +
  scale_color_gradient2(low = "#2b83ba", mid = "#FFFFB3", high = "#d7191c",
                        midpoint = 1, limits = c(0, 2), oob = scales::squish) +
  geom_hline(yintercept = -log10(cutoff_p), color = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "solid") +
  labs(title = 'The Replicated HSV peptides for Prevalent Diseases',
       y = expression('-log'[10]*'('*italic(P)*'-value)'),
       color = "OR") +
  theme_bw() +
  theme(
    panel.border  = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text  = element_text(size = 8),
    plot.title   = element_text(hjust = 0.5, size = 11),
    strip.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.text.x  = element_text(angle = 70, hjust = 1, size = 8),
    axis.text.y  = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  )

p.all <- p.all + ggrepel::geom_label_repel(
  data = subset(res_wide, P_all < 1e-04),
  aes(label = taxon_species),
  size = 2, segment.size = 0.1, direction = "both",
  segment.color = 'black', max.overlaps = 20
)

ggsave(file.path(opt$out_dir, "figure5_or_pre_all_virus.png"), p.all, width = 12, height = 7, dpi = 300, bg='white')

# OR significance comparison
p.or <- ggplot(res_wide, aes(x = disease, y = -log10(P_all))) +
  facet_grid(. ~ Category, scales = "free_x") +
  geom_point(aes(color = Sig_Status), size = 2) +
  scale_color_manual(values = c("Both Significant" = "red", 
                                "MGBB-LEC Only" = "blue", 
                                "MGBB-LLF Only" = "green", 
                                "Not Significant" = "grey80")) +
  geom_hline(yintercept = -log10(cutoff_p), color = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "solid") +
  labs(title = 'The Replicated HSV peptides for Prevalent Diseases',
       y = expression('-log'[10]*'('*italic(P)*'-value)'),
       color = "OR") +
  theme_bw() +
  theme(
    panel.border  = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text  = element_text(size = 8),
    plot.title   = element_text(hjust = 0.5, size = 11),
    strip.text.x = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.text.x  = element_text(angle = 70, hjust = 1, size = 8),
    axis.text.y  = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  )

p.or <- p.or + ggrepel::geom_label_repel(
  data = subset(res_wide, P_all < 1e-04),
  aes(label = taxon_species),
  size = 2, segment.size = 0.1, direction = "both",
  segment.color = 'black', max.overlaps = 20
)
ggsave(file.path(opt$out_dir, "figure5_or_comparison_glm_all_virus.png"), p.or, width = 12, height = 7, dpi = 300, bg='white')

# Beta comparison Plot
p.beta = ggplot(res_wide, aes(x = Beta_lec, y = Beta_llf, color = Sig_Status)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted") + # Perfect replication line
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Both Significant" = "red", 
                                "MGBB-LEC Only" = "blue", 
                                "MGBB-LLF Only" = "green", 
                                "Not Significant" = "grey80")) +
  theme_minimal() +
  labs(title = "Cohort comparison between MGBB-LLF and MGBB-LEC",
       subtitle = "Checking for consistency of effect direction",
       x = "Effect Size (MGBB-LEC)", 
       y = "Effect Size (MGBB-LLF)") +
  theme(
    panel.border  = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text  = element_text(size = 10),
    plot.title   = element_text(size = 12),
    axis.title.x = element_text(size = 10),
    axis.text.x  = element_text(size = 8),
    axis.text.y  = element_text(size = 10),
    axis.title.y = element_text(size = 10)
  ) +
  # Optional: Label only the top significant hits to avoid clutter
  ggrepel::geom_label_repel(data = res_wide[Sig_Status == "Both Significant"], 
                            aes(label = paste(disease, taxon_species, sep="-")),
                            size = 3, segment.size = 0.1, direction = "both",
                            segment.color = 'black', max.overlaps = 20)

ggsave(file.path(opt$out_dir, "figure5_beta_glm_all_virus.png"), p.beta, width = 12, height = 10, dpi = 300, bg='white')

