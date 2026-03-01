#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse); library(data.table); library(dplyr); library(glmnet)
  library(pROC); library(ggplot2); library(ggrepel); library(scales)
  library(corrplot); library(purrr); library(stringr); library(tidyr);library(foreach);library(doParallel)
  source("R/utils_common.R")
  source("R/utils_qtl.R")
})

# ----------------------------- CLI -----------------------------
option_list <- list(
  make_option("--train_vi", type="character", help="hsv_promax_bin_MGBB-LLF.tsv"),
  make_option("--train_hu", type="character", help="human_fl_bin_MGBB-LLF.tsv"),
  make_option("--train_vi_anno", type="character", help="IB1007_VirSIGHT_Promax_Hits_Fold-Over-Background.csv"),
  make_option("--train_hu_anno", type="character", help="IB1007_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv"),
  make_option("--test_vi", type="character", help="hsv_promax_bin_MGBB-LEC.tsv"),
  make_option("--test_hu", type="character", help="human_fl_bin_MGBB-LEC.tsv"),
  make_option("--test_hu", type="character", help="hsv_bin_fchange_rep_llf.tsv"),
  make_option("--out_dir", type="character", default="results/Figure4", help="Output dir")
)

opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$out_dir, showWarnings=FALSE, recursive=TRUE)

# ----------------------------- IO ------------------------------
train_vi <- fread(opt$train_vi); train_hu <- fread(opt$train_hu)
test_vi  <- fread(opt$test_vi);  test_hu  <- fread(opt$test_hu)

rep_hsv <- fread(opt$rep_hsv)

train_hu_anno <- fread(opt$train_hu_anno)
train_vi_anno <- fread(opt$train_vi_anno)

if (!"var_id" %in% names(train_hu_anno)) train_hu_anno[, var_id := paste0("h", seq_len(.N))]
train_vi_anno$taxon_species <- species_label_map(train_vi_anno$taxon_species)

# ------------------ Focus antibodies (Figure 4 core) -----------
rep_ab_tbl = unique(rep_hsv[,c('antibody','gene_symbol')])

# Shared-subject join
merge_train <- train_hu %>% right_join(train_vi, by="Subject_Id") %>% as.data.table()
merge_test  <- test_hu  %>% right_join(test_vi,  by="Subject_Id") %>% as.data.table()

# virus species sets per short label
species_sets <- lapply(species_list, function(sp) {
  acc <- train_vi_anno %>% filter(taxon_species == sp) %>% pull(UniProt_acc) %>%
    intersect(colnames(train_vi)) %>% intersect(colnames(test_vi))
  tibble::tibble(species = sp, acc = list(acc))
}) %>% bind_rows()

# 1. Setup Parallel Backend
# Leave one core free for OS responsiveness
num_cores <- parallel::detectCores() - 1 
registerDoParallel(cores = num_cores)

print(paste("Running on", num_cores, "cores"))

# 2. Run Parallel Loop over Species (Outer Loop)
# We parallelize the outer loop so each core handles one species and its matrix generation
results_list.lec <- foreach(i = seq_len(nrow(species_sets))) %dopar% {
  
  sp  <- species_sets$species[i]
  acc <- species_sets$acc[[i]]

  if (length(acc) == 0) return(NULL)
  
  # --- OPTIMIZATION A: Pre-compute Matrix ONCE per species ---
  valid_cols <- intersect(acc, colnames(merge_train))
  if (length(valid_cols) == 0) return(NULL)
  
  # Create matrices (Using -1 to remove intercept if you prefer glmnet to handle it, 
  # though keeping default is fine if consistent)
  train_x <- model.matrix(~ ., merge_train[, ..valid_cols])
  test_x  <- model.matrix(~ ., merge_test[, ..valid_cols])
  
  # Local list to store results for this species
  species_results <- list()
  
  # Iterate antibodies (Inner Loop - runs sequentially on each core)
  for (j in seq_len(nrow(rep_ab_tbl))) {
    ab   <- rep_ab_tbl$antibody[j]
    gene <- rep_ab_tbl$gene_symbol[j]
    
    # Skip if antibody missing
    if (!(ab %in% colnames(merge_train)) || !(ab %in% colnames(merge_test))) next
    
    train_y <- merge_train[[ab]]
    test_y  <- merge_test[[ab]]
    
    # Skip if target has no variance (all 0s or all 1s)
    if (length(unique(train_y)) < 2) next
    
    # --- OPTIMIZATION B: Faster CV ---
    set.seed(123)
    cv_fit <- tryCatch(
      cv.glmnet(train_x, train_y, family="binomial", alpha=1, nfolds=10),
      error = function(e) NULL
    )
    
    if (is.null(cv_fit)) next
    
    # Prediction
    pred_prob <- as.numeric(predict(cv_fit, newx=test_x, s="lambda.min", type="response"))
    
    # ROC Calculation
    roc_obj <- tryCatch(roc(test_y, pred_prob, quiet = TRUE), error=function(e) NULL)
    
    if (!is.null(roc_obj)) {
      auc_val <- as.numeric(auc(roc_obj))
      
      # We create the summary row immediately to save memory
      res_row <- data.frame(
        Antibody = ab,
        Gene = gene,
        Virus = sp,
        AUC = auc_val,
        Lambda = cv_fit$lambda.min
      )
      
      roc_coords <- coords(roc_obj, x="all", ret=c("specificity","sensitivity")) %>%
        as.data.frame() %>% 
        mutate(FPR = 1 - specificity, TPR = sensitivity, Antibody=ab, Virus=sp)
      
      species_results[[paste0(ab, "_", sp)]] <- list(summary = res_row, roc = roc_coords)
    }
  }
  
  # Return the list of results for this species
  return(species_results)
}

# 3. Post-process: Unpack the parallel results
# Flatten the list of lists
flat_results.lec <- unlist(results_list.lec, recursive = FALSE)

# Extract DataFrames
if (!is.null(flat_results.lec)) {
  auc_tab_df.lec <- bind_rows(lapply(flat_results.lec, `[[`, "summary"))
  roc_all_df.lec <- bind_rows(lapply(flat_results.lec, `[[`, "roc"))
}

roc_all_df.lec = left_join(roc_all_df.lec, auc_tab_df.lec, by=c("Antibody", 'Virus'))

auc_tab_df.lec1 <- auc_tab_df.lec %>%  
  group_by(Gene, Virus) %>%
  slice_max(order_by = AUC, with_ties = FALSE) %>%
  ungroup()
auc_tab_df.lec1 = left_join(auc_tab_df.lec1, train_hu_anno[,c('var_id','UniProt_acc')], by=c('Antibody'='var_id'))
auc_tab_df.lec1$Gene <- gsub(" \\s*\\([^\\)]+\\)", "",  auc_tab_df.lec1$Gene)
auc_tab_df.lec1$Gene <- gsub("GPRASP3,ARMCX5-GPRASP2", "GPRASP3", auc_tab_df.lec1$Gene)

saveRDS(flat_results.lec, file = file.path(opt$out_dir, "lec_test_lasso.rds"))
fwrite(roc_all_df.lec, file.path(opt$out_dir, "lec_test_roc_summary.csv"))
fwrite(auc_tab_df.lec1, file.path(opt$out_dir, "lec_test_auc_summary.csv"))
# flat_results.lec = readRDS(file = file.path(opt$out_dir, "lec_test_lasso.rds"))

ab_tbl = subset(auc_tab_df.lec1[,c('Antibody', 'Gene', 'Virus','AUC')], auc_tab_df.lec1$AUC > 0.85) #17
fwrite(ab_tbl, file.path(opt$out_dir, "lec_test_auc_summary_sig.csv"))

# ------------------ LLF only ------------------------------------------------

# 1. Parallel Loop over Species in LLF only
results_list.llf <- foreach(i = seq_len(nrow(species_sets))) %dopar% {
  
  sp  <- species_sets$species[i]
  acc <- species_sets$acc[[i]]
  print(sp)
  
  # Basic checks
  if (length(acc) == 0) return(NULL)
  valid_cols <- intersect(acc, colnames(merge_train)) # Assuming 'merge_train' is total data now
  if (length(valid_cols) == 0) return(NULL)
  
  # --- OPTIMIZATION 1: Create ONE Design Matrix for all antibodies of this species ---
  # We build the matrix on the TOTAL data first
  full_x <- model.matrix(~ ., merge_train[, ..valid_cols])
  
  # --- OPTIMIZATION 2: Perform 80:20 Split ---
  # Create indices: 80% Train, 20% Test
  n_samples <- nrow(merge_train)
  # Set seed based on iteration to ensure reproducible splits in parallel
  set.seed(123)
  train_idx <- sample(seq_len(n_samples), size = floor(0.8 * n_samples))
  
  # Split the matrix immediately (so we don't do it inside the antibody loop)
  train_x_mat <- full_x[train_idx, ]
  test_x_mat  <- full_x[-train_idx, ]
  
  species_results <- list()
  
  # --- Inner Loop: Iterate Antibodies ---
  for (j in seq_len(nrow(rep_ab_tbl))) {
    ab   <- rep_ab_tbl$antibody[j]
    gene <- rep_ab_tbl$gene_symbol[j]

    # Check if antibody exists
    if (!(ab %in% colnames(merge_train))) next
    
    # Get the response vector (y) for all samples
    full_y <- merge_train[[ab]]
    
    # Split y using the SAME indices as the matrix
    train_y_vec <- full_y[train_idx]
    test_y_vec  <- full_y[-train_idx]
    
    # --- Safety Checks ---
    # 1. Skip if training data has no variation (all 0s or all 1s)
    if (length(unique(train_y_vec)) < 2) next
    # 2. Skip if test set has no positives (cannot calculate ROC)
    if (sum(test_y_vec == 1) == 0 || sum(test_y_vec == 0) == 0) next
    
    # --- Lasso with Faster CV ---
    set.seed(123)
    cv_fit <- tryCatch(
      cv.glmnet(train_x_mat, train_y_vec, family="binomial", alpha=1, nfolds=10),
      error = function(e) NULL
    )
    
    if (is.null(cv_fit)) next
    
    # Prediction on Test Set
    pred_prob <- as.numeric(predict(cv_fit, newx=test_x_mat, s="lambda.min", type="response"))
    
    # ROC / AUC
    roc_obj <- tryCatch(roc(test_y_vec, pred_prob, quiet = TRUE), error=function(e) NULL)
    
    if (!is.null(roc_obj)) {
      auc_val <- as.numeric(auc(roc_obj))
      
      # Save summary row
      res_row <- data.frame(
        Antibody = ab,
        Gene = gene,
        Virus = sp,
        AUC = auc_val,
        Lambda = cv_fit$lambda.min
      )
      
      # Save ROC coordinates (Optional - remove if file size is too big)
      roc_coords <- coords(roc_obj, x="all", ret=c("specificity","sensitivity")) %>%
        as.data.frame() %>% 
        mutate(FPR = 1 - specificity, TPR = sensitivity, Antibody=ab, Virus=sp)
      
      species_results[[paste0(ab, "_", sp)]] <- list(summary = res_row, roc = roc_coords)
    }
  }
  
  return(species_results)
}

# 3. Clean up Parallel Backend
stopImplicitCluster()

# 4. Consolidate Results
flat_results.llf <- unlist(results_list.llf, recursive = FALSE)

auc_tab_df.llf <- bind_rows(lapply(flat_results.llf, `[[`, "summary"))
roc_all_df.llf <- bind_rows(lapply(flat_results.llf, `[[`, "roc"))

roc_all_df.llf = left_join(roc_all_df.llf, auc_tab_df.llf, by=c("Antibody", 'Virus'))

auc_tab_df.llf1 <- auc_tab_df.llf %>%  
  group_by(Gene, Virus) %>%
  slice_max(order_by = AUC, with_ties = FALSE) %>%
  ungroup()
auc_tab_df.llf1 = left_join(auc_tab_df.llf1, train_hu_anno[,c('var_id','UniProt_acc')], by=c('Antibody'='var_id'))
auc_tab_df.llf1$Gene <- gsub(" \\s*\\([^\\)]+\\)", "",  auc_tab_df.llf1$Gene)
auc_tab_df.llf1$Gene <- gsub("GPRASP3,ARMCX5-GPRASP2", "GPRASP3", auc_tab_df.llf1$Gene)

saveRDS(flat_results.llf, file = file.path(opt$out_dir, "llf_test_lasso.rds"))
# flat_results.llf = readRDS(file = file.path(opt$out_dir, "llf_test_lasso.rds"))

fwrite(roc_all_df.llf, file.path(opt$out_dir, "llf_test_roc_summary.csv"))
fwrite(auc_tab_df.llf, file.path(opt$out_dir, "llf_test_auc_summary.csv"))

ab_tbl.llf = subset(auc_tab_df.llf[,c('Antibody', 'Gene', 'Virus','AUC')], auc_tab_df.llf$AUC > 0.80) #28
fwrite(ab_tbl.llf, file.path(opt$out_dir, "llf_test_auc_summary_sig.csv"))

auc_merge = inner_join(auc_tab_df.lec1, auc_tab_df.llf, by=c('Antibody', 'Gene', 'Virus'))

fwrite(auc_merge, file.path(opt$out_dir, "auc_summary_merge.csv"))
##################################################################################
# AUC >0.85 Gene-HSV links

roc_all_df.lec$Virus = factor(roc_all_df.lec$Virus, levels = c("HSV-1","HSV-2","VZV","CMV","HHV-6A","HHV-6B",'HHV-7',"EBV",'HHV-8'))

if (nrow(roc_all_df.lec)) {
  for (ab in unique(ab_tbl$Antibody)) {
    sub  <- roc_all_df.lec %>% filter(Antibody==ab)
    gene <- unique(sub$Gene)
    p <- ggplot(sub, aes(FPR, TPR, color = Virus)) +
      geom_line(size=1) + geom_abline(linetype="dashed") +
      labs(title = paste0("", gene, " MaxAUC=", sprintf("%.3f", max(sub$AUC))),
           x="False Positive Rate", y="True Positive Rate") +
      theme_minimal(base_size=12)
    ggsave(file.path(opt$out_dir, paste0("roc_", gene, "_lec.png")), p, width=5.5, height=5, dpi=300, bg = "white")
  }
}

roc_all_df.llf$Virus = factor(roc_all_df.llf$Virus, levels = c("HSV-1","HSV-2","VZV","CMV","HHV-6A","HHV-6B",'HHV-7',"EBV",'HHV-8'))

if (nrow(roc_all_df.llf)) {
  for (ab in unique(ab_tbl$Antibody)) {
    sub  <- roc_all_df.llf %>% filter(Antibody==ab)
    gene <- unique(sub$Gene)
    p <- ggplot(sub, aes(FPR, TPR, color = Virus)) +
      geom_line(size=1) + geom_abline(linetype="dashed") +
      labs(title = paste0("", gene, " MaxAUC=", sprintf("%.3f", max(sub$AUC))),
           x="False Positive Rate", y="True Positive Rate") +
      theme_minimal(base_size=12)
    ggsave(file.path(opt$out_dir, paste0("roc_", gene, "_llf.png")), p, width=5.5, height=5, dpi=300, bg = "white")
  }
}
# Supplementary Figure 1 across HSVs
# Top antibodies per virus

for (ii in c(1:9)) {
  virus_name = species_list[ii]
  paint = color[ii]
  message("Processing: ", virus_name)
  
  subgroup <- subset(auc_tab_df.lec, Virus == virus_name)
  subgroup1 <- subgroup %>%
    group_by(Gene) %>%
    slice_max(order_by = AUC, with_ties = FALSE) %>%
    ungroup()
  
  subgroup1$Gene <- gsub(" \\s*\\([^\\)]+\\)", "",  subgroup1$Gene)
  subgroup1 <- subgroup1[order(subgroup1$AUC, decreasing = F),]
  subgroup1$Gene <- factor(subgroup1$Gene, levels = subgroup1$Gene)

  p <- ggplot(subgroup1, aes(x = Gene, y = AUC)) +
    geom_bar(stat = "identity", fill=paint) +
    coord_flip() +
    ylim(0, 1) +  
    theme_minimal(base_size = 12) +
    labs(
      title = paste(virus_name),
      x = "Gene Symbol",
      y = "AUC"
    ) +
    theme(plot.title = element_text(face = "bold", size = 14), axis.title  = element_text(size = 12), axis.text   = element_text(size = 11))
  
  print(p)
  
  ggsave(file.path(opt$out_dir, paste0("auc_by_",virus_name,"_lec.png")), plot = p, width = 9, height = 14, dpi = 300, bg = "white")
}

for (ii in c(1:9)) {
  virus_name = species_list[ii]
  paint = color[ii]
  message("Processing: ", virus_name)
  
  subgroup <- subset(auc_tab_df.llf, Virus == virus_name)
  subgroup1 <- subgroup %>%
    group_by(Gene) %>%
    slice_max(order_by = AUC, with_ties = FALSE) %>%
    ungroup()
  
  subgroup1$Gene <- gsub(" \\s*\\([^\\)]+\\)", "",  subgroup1$Gene)
  subgroup1$Gene <- gsub("GPRASP3,ARMCX5-GPRASP2", "GPRASP3", subgroup1$Gene)
  subgroup1 <- subgroup1[order(subgroup1$AUC, decreasing = F),]
  subgroup1$Gene <- factor(subgroup1$Gene, levels = subgroup1$Gene)
  
  p <- ggplot(subgroup1, aes(x = Gene, y = AUC)) +
    geom_bar(stat = "identity", fill=paint) +
    coord_flip() +
    ylim(0, 1) +  
    theme_minimal(base_size = 12) +
    labs(
      title = paste(virus_name),
      x = "Gene Symbol",
      y = "AUC"
    ) +
    theme(plot.title = element_text(face = "bold", size = 14), axis.title  = element_text(size = 12), axis.text   = element_text(size = 11))
  
  print(p)
  
  ggsave(file.path(opt$out_dir, paste0("auc_by_",virus_name,"_llf.png")), plot = p, width = 9, height = 14, dpi = 300, bg = "white")
}

# Correlation among key Abs (train set)
ab_tbl = subset(auc_merge, auc_merge$AUC.x>0.85 & auc_merge$AUC.y>0.8)
ab_tbl = ab_tbl[order(ab_tbl$Virus, ab_tbl$Gene),]
keep_abs <- as.factor(ab_tbl$Antibody)

if (length(keep_abs) >= 2) {
  cor_mat <- cor(merge_train[, ..keep_abs], method="pearson", use="pairwise.complete.obs")
  colnames(cor_mat) <- ab_tbl$Gene[match(colnames(cor_mat), ab_tbl$Antibody)]
  rownames(cor_mat) <- colnames(cor_mat)
  ragg::agg_png(
    file.path(opt$out_dir, "figure4F_correlation_llf.png"),
    width = 2100, height = 2100, res = 300
  )
  corrplot(cor_mat, method="color", type="upper", title = "", mar = c(0, 0, 2, 0),
           number.cex=0.7, tl.cex=1.0, tl.col="black", tl.srt=45, addCoef.col="black")
  dev.off()
}



# ---- main function (drop-in) ----
run_lasso_feature_analysis <- function(
    virus_name,              # "HSV-1" or long name
    antibody_id,             # e.g., "h6808"
    gene_id,                 # e.g., "PHLDA1" (for titles)
    train_vi_anno,           # annotation table with at least UniProt_acc, taxon_species, product
    train_vi, train_hu,      # MGBB/LLF binary matrices with Sample_Id and features
    test_vi, test_hu,        # LEC binary matrices with Sample_Id and features
    out_dir,                 # output directory
    color_map
) {
  message("\nRunning detailed feature LASSO for ", virus_name, " / ", antibody_id, " / ", gene_id)
  
  # standardize species labels in the annotation
  train_vi_anno <- train_vi_anno %>% mutate(taxon_species = species_label_map(taxon_species))

  # 1) collect UniProt IDs for this virus present in BOTH train and test
  virus_pro <- train_vi_anno %>%
    filter(taxon_species == virus_name) %>%
    pull(UniProt_acc) %>%
    intersect(colnames(train_vi)) %>%
    intersect(colnames(test_vi))
  
  if (!length(virus_pro)) {
    message("No viral proteins found for ", virus_name, " that are shared in train and test.")
    return(NULL)
  }
  
  # 2) build train & test joined frames
  train_all <- train_hu %>% right_join(train_vi, by = "Subject_Id")
  test_all  <- test_hu  %>% right_join(test_vi,  by = "Subject_Id")
  
  needed_cols <- c(virus_pro, antibody_id)
  if (!all(needed_cols %in% colnames(train_all))) {
    message("Missing columns in training data for ", virus_name, " / ", antibody_id)
    return(NULL)
  }
  if (!all(needed_cols %in% colnames(test_all))) {
    message("Missing columns in test data for ", virus_name, " / ", antibody_id)
    return(NULL)
  }
  
  # keep only required columns
  df_train <- train_all %>% select(all_of(needed_cols))
  df_test  <- test_all  %>% select(all_of(needed_cols))
  y_tr <- df_train[[antibody_id]]
  y_te <- df_test[[antibody_id]]
  
  # 3) model matrices
  # combine to build a consistent design; then split back
  df_comb <- rbind(df_train, df_test)
  X <- model.matrix(
    reformulate(termlabels = setdiff(colnames(df_comb), antibody_id), response = antibody_id),
    data = df_comb
  )
  y <- df_comb[[antibody_id]]
  
  n_tr <- nrow(df_train)
  X_train <- X[seq_len(n_tr), , drop = FALSE]
  y_train <- y[seq_len(n_tr)]
  X_test  <- X[-seq_len(n_tr), , drop = FALSE]
  y_test  <- y[-seq_len(n_tr)]
  
  # glmnet can handle the intercept internally
  if ("(Intercept)" %in% colnames(X_train)) {
    Xi <- which(colnames(X_train) == "(Intercept)")
    if (length(Xi) == 1) {
      X_train <- X_train[, -Xi, drop = FALSE]
      X_test  <- X_test[,  -Xi, drop = FALSE]
    }
  }
  
  # 4) LASSO fit
  set.seed(123)
  cv_fit <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 1, nfolds = 10)
  best_lambda <- cv_fit$lambda.min
  
  # non-zero features at lambda.min
  coef_mat <- as.matrix(coef(cv_fit, s = best_lambda))
  nz <- which(as.numeric(coef_mat) != 0)
  vars <- rownames(coef_mat)[nz]
  betas <- as.numeric(coef_mat[nz, 1])
  # remove intercept if present
  keep <- vars != "(Intercept)"
  vars <- vars[keep]; betas <- betas[keep]
  
  keys <- tibble::tibble(variables = vars, beta = betas)
  
  # 5) annotate with product names
  keys1 <- keys %>%
    inner_join(train_vi_anno[, c("UniProt_acc","product")], by = c("variables" = "UniProt_acc")) %>%
    mutate(product = as.factor(product))
  
  # 6) OLS p-values (post-selection; treat cautiously)
  p_df <- NULL
  if (length(vars) >= 1) {
    Xsel <- as.data.frame(X_train[, vars, drop = FALSE])
    ols <- tryCatch(summary(lm(y_train ~ ., data = Xsel)), error = function(e) NULL)
    if (!is.null(ols)) {
      p_df <- data.frame(
        variables = rownames(ols$coefficients),
        P_Value   = ols$coefficients[,4],
        row.names = NULL,
        check.names = FALSE
      )
      # drop intercept if it shows up
      p_df <- p_df[p_df$variables != "(Intercept)", , drop = FALSE]
    }
  }
  keys2 <- if (!is.null(p_df) && nrow(p_df)) {
    left_join(keys1, p_df, by = "variables")
  } else {
    mutate(keys1, P_Value = NA_real_)
  }
  keys2 <- keys2 %>%
    mutate(sig_marker = dplyr::case_when(
      P_Value < 0.001 ~ "***",
      P_Value < 0.01  ~ "**",
      P_Value < 0.05  ~ "*",
      TRUE ~ ""
    ))
  
  keys3 = subset(keys2, keys2$P_Value < 0.05)
  
  for (pat in names(clean_replace)) {
    keys3$product <- gsub(pat, clean_replace[[pat]], keys3$product, fixed = TRUE)
  }
  
  # 7) LEC per-variable 2×2 metrics (sensitivity/specificity/AUC_1threshold)
  #    Evaluate *each* selected viral peptide individually on LEC
  results_lec <- lapply(keys3$variables, function(v) {
    if (!v %in% colnames(test_vi)) return(NULL)
    # make sure they're 0/1 numeric
    a <- as.integer(test_hu[[antibody_id]])
    b <- as.integer(test_vi[[v]])
    if (any(is.na(a)) || any(is.na(b))) {
      ok <- which(!is.na(a) & !is.na(b))
      a <- a[ok]; b <- b[ok]
    }
    tt <- table(a, b)
    if (!identical(dim(tt), c(2L, 2L))) return(NULL)
    TN <- as.numeric(tt[1,1]); FP <- as.numeric(tt[1,2])
    FN <- as.numeric(tt[2,1]); TP <- as.numeric(tt[2,2])
    sens <- TP/(TP+FN); spec <- TN/(TN+FP)
    auc1 <- (sens + spec)/2
    pv <- tryCatch(chisq.test(tt)$p.value, error = function(e) NA_real_)
    data.frame(Variable=v, Sensitivity=sens, Specificity=spec,
               AUC_1threshold=auc1, P.value=pv, hid=antibody_id, Gene=gene_id, stringsAsFactors = FALSE)
  })
  results_lec <- data.table::rbindlist(results_lec, fill = TRUE, use.names = TRUE)
  results_lec <- results_lec[order(-AUC_1threshold), ]
  results_lec1 <- left_join(results_lec, train_vi_anno[, c("UniProt_acc","product", "taxon_species")],
                            by = c("Variable"="UniProt_acc"))
  
  # save the table
  data.table::fwrite(results_lec1,
                     file.path(out_dir, paste0("figure4_", virus_name, "_", gene_id, "_LEC_metrics.csv")),
                     sep = ",")
  
  # 9) LEC scatter (Sensitivity vs Specificity, colored by AUC_1threshold)
  plot_data <- results_lec %>% filter(!is.na(Sensitivity), !is.na(Specificity))
  p_roc <- ggplot(plot_data, aes(x = Sensitivity, y = Specificity,
                                 color = AUC_1threshold, label = Variable)) +
    geom_point(size = 2.8) +
    ggrepel::geom_text_repel(size = 2.6, show.legend = FALSE, max.overlaps = 40) +
    scale_color_viridis_c(option = "plasma", limits = c(0.5, 1), name = expression(AUC)) +
    theme_minimal() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "Sensitivity", y = "Specificity", title = paste0("",virus_name," peptides for ",gene_id,"")) +
    theme(legend.text = element_text(size = 8),
          plot.title  = element_text(hjust = 0.5, size = 12),
          axis.title  = element_text(size = 10),
          axis.text   = element_text(size = 8))
  
  ragg::agg_png(file.path(out_dir, paste0("figure4_", virus_name, "_", gene_id, "_ROC.png")),
                width = 1000, height = 1000, res = 300)
  print(p_roc); dev.off()
  
  invisible(list(coef = keys2, lec_metrics = results_lec1))
}

# Run for the main 5 panels
pairs_to_run <- split(ab_tbl, seq(nrow(ab_tbl)))

for (pair in pairs_to_run) {
  run_lasso_feature_analysis(pair$Virus, pair$Antibody, pair$Gene,
                             train_vi_anno = train_vi_anno,
                             train_vi = train_vi, train_hu = train_hu,
                             test_vi = test_vi,   test_hu  = test_hu,
                             out_dir = opt$out_dir,
                             color_map = color_map)
}

system(paste0("awk '(NR == 1) || (FNR > 1)' ",opt$out_dir,"figure4*_metrics.csv  > ",opt$out_dir,"auc_virus_metrics.csv"))
