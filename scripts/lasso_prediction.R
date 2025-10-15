#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(dplyr); library(glmnet)
  library(pROC); library(ggplot2); library(ggrepel); library(scales)
  library(corrplot); library(purrr); library(stringr); library(tidyr)
  source("R/utils_common.R")
  source("R/utils_utils.R")
})

# ----------------------------- CLI -----------------------------
option_list <- list(
  make_option("--train_vi", type="character", help="MGBB LLF virus_proteins_binary.tsv"),
  make_option("--train_hu", type="character", help="MGBB LLF human_fl_binary.tsv"),
  make_option("--train_vi_anno", type="character", help="LLF VirSIGHT Promax FoB CSV"),
  make_option("--train_hu_anno", type="character", help="LLF HuSIGHT FullLength FoB CSV"),
  make_option("--test_vi", type="character", help="ABC virus binary matrix"),
  make_option("--test_hu", type="character", help="ABC human binary matrix"),
  # Figure 4 focus set (specific antibodies to showcase ROC, coefs, corr)
  make_option("--ab_list", type="character",
              help="Comma sep gene:hid pairs, e.g. 'PHLDA1:h6808,ZNF550:h6329,IQCB1:h3143,DNAJC12:h8069,P3H4:h8250'"),
  make_option("--species_map", type="character",
              help="Semicolon sep 'Human alphaherpesvirus 1=HSV-1;Human betaherpesvirus 5=CMV;Human gammaherpesvirus 4=EBV'"),
  make_option("--out_dir", type="character", default="results/Figure4", help="Output dir"),
  make_option("--make_combined_panel", action="store_true", default=FALSE,
              help="If set, combine 9 HSV barplots into a 3x3 panel PNG")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$out_dir, showWarnings=FALSE, recursive=TRUE)

# ----------------------------- IO ------------------------------
train_vi <- fread(opt$train_vi); train_hu <- fread(opt$train_hu)
test_vi  <- fread(opt$test_vi);  test_hu  <- fread(opt$test_hu)

rep_hsv <- fread(opt$rep_hsv)

train_hu_anno <- fread(opt$train_hu_anno)
train_vi_anno <- fread(opt$train_vi_anno)

test_hu_anno <- fread(opt$test_hu_anno)
test_vi_anno <- fread(opt$test_vi_anno)

if (!"var_id" %in% names(train_hu_anno)) train_hu_anno[, var_id := paste0("h", seq_len(.N))]

bin_from_fob <- function(x) as.integer(!is.na(as.numeric(x)) & as.numeric(x) > 1)
test_vi_anno$taxon_species <- species_label_map(test_vi_anno$taxon_species)
test_vi_anno1 <- test_vi_anno[taxon_species %in% species_list]

sample_cols_vi <- grep('^10', names(test_vi_anno1), value = TRUE)
virus_peptide <- grep('Sample_ID', names(train_vi), invert = TRUE, value = TRUE)

test_vi_anno2 <- as.data.frame(test_vi_anno1[, ..sample_cols_vi])
test_vi_anno2[] <- lapply(test_vi_anno2, bin_from_fob)
test_vi_anno2 <- as.data.frame(t(test_vi_anno2))
colnames(test_vi_anno2) <- test_vi_anno1$UniProt_acc
rownames(test_vi_anno2) <- NULL
test_vi = test_vi_anno2[,virus_peptide]
test_vi$Sample_ID = sample_cols_vi
test_hu$Sample_ID = as.character(test_hu$Sample_ID)

# ------------------ Focus antibodies (Figure 4 core) -----------
ab_pairs <- str_split(opt$ab_list, ",")[[1]]
ab_tbl   <- tibble::tibble(
  gene = sub(":.*$", "", ab_pairs),
  hid  = paste0(sub("^.*:", "", ab_pairs))
)

rep_ab_tbl = unique(rep_hsv[,c('antibody','gene_symbol')])

# Shared-subject join
merge_train <- train_hu %>% right_join(train_vi, by="Sample_ID") %>% as.data.table()
merge_test  <- test_hu  %>% right_join(test_vi,  by="Sample_ID") %>% as.data.table()

# virus species sets per short label
species_sets <- lapply(species_list, function(sp) {
  acc <- train_vi_anno %>% filter(taxon_species == sp) %>% pull(UniProt_acc) %>%
    intersect(colnames(train_vi)) %>% intersect(colnames(test_vi))
  tibble::tibble(species = sp, acc = list(acc))
}) %>% bind_rows()

fit_lasso_eval <- function(species_acc, ab_id, df_train, df_test, out_prefix) {
  acc <- species_acc
  if (length(acc) == 0) return(NULL)
  if (!(ab_id %in% colnames(df_train)) || !(ab_id %in% colnames(df_test))) return(NULL)

  train_x <- model.matrix(~ ., df_train[, ..acc])
  train_y <- df_train[[ab_id]]
  test_x  <- model.matrix(~ ., df_test[,  ..acc])
  test_y  <- df_test[[ab_id]]
  
  set.seed(123)
  cv_fit   <- cv.glmnet(train_x, train_y, family="binomial", alpha=1)
  pred_prob<- as.numeric(predict(cv_fit, newx=test_x, s=cv_fit$lambda.min, type="response"))
  
  roc_obj <- tryCatch(roc(test_y, pred_prob, quiet = TRUE), error=function(e) NULL)
  if (is.null(roc_obj)) return(NULL)
  
  auc_val <- as.numeric(auc(roc_obj))
  roc_df  <- coords(roc_obj, x="all", ret=c("specificity","sensitivity")) %>%
    as.data.frame() %>% mutate(FPR = 1 - specificity, TPR = sensitivity)
  
  list(roc = roc_df, auc = auc_val, lambda = cv_fit$lambda.min, cv = cv_fit)
}

# ----------------------- Figure 4: ROC etc. --------------------
roc_all <- list(); auc_tab <- list()

for (i in seq_len(nrow(species_sets))) {
  sp  <- species_sets$species[i]
  acc <- species_sets$acc[[i]]
  for (j in seq_len(nrow(rep_ab_tbl))) {
    ab   <- rep_ab_tbl$antibody[j]
    gene <- rep_ab_tbl$gene_symbol[j]
    res <- fit_lasso_eval(acc, ab, merge_train, merge_test, paste0(sp,"_",ab))
    if (is.null(res)) next
    roc_all[[paste0(ab,"_",sp)]] <- res$roc %>% mutate(Antibody=ab, Gene=gene, Virus=sp, AUC=res$auc)
    auc_tab[[paste0(ab,"_",sp)]] <- data.frame(Antibody=ab, Gene=gene, Virus=sp, AUC=res$auc)
  }
}

roc_all_df <- bind_rows(roc_all)
auc_tab_df <- bind_rows(auc_tab)

fwrite(roc_all_df, file.path(opt$out_dir, "roc_summary_ABCtest.csv"))
fwrite(auc_tab_df, file.path(opt$out_dir, "auc_summary_ABCtest.csv"))


roc_all_df$Virus = factor(roc_all_df$Virus, levels = c("HSV-1","HSV-2","VZV","CMV","HHV-6A","HHV-6B",'HHV-7',"EBV",'HHV-8'))

if (nrow(roc_all_df)) {
  for (ab in unique(ab_tbl$hid)) {
    sub  <- roc_all_df %>% filter(Antibody==ab)
    gene <- unique(sub$Gene)
    p <- ggplot(sub, aes(FPR, TPR, color = Virus)) +
      geom_line(size=1) + geom_abline(linetype="dashed") +
      labs(title = paste0("ROC for ", gene, " (ABC test)  MaxAUC=", sprintf("%.3f", max(sub$AUC))),
           x="False Positive Rate", y="True Positive Rate") +
      theme_minimal(base_size=12)
    ggsave(file.path(opt$out_dir, paste0("roc_", ab, ".png")), p, width=5.5, height=5, dpi=300)
  }
}


# Supplementary Figure 1 across HSVs
for (ii in c(1:9)) {
  virus_name = species_list[ii]
  paint = color[ii]
  message("Processing: ", virus_name)
  
  subgroup <- subset(auc_tab_df, Virus == virus_name)
  subgroup$Gene <- make.unique(as.character(subgroup$Gene))
  subgroup <- subgroup[order(subgroup$AUC, decreasing = F),]
  subgroup$Gene <- factor(subgroup$Gene, levels = subgroup$Gene)
  
  p <- ggplot(subgroup, aes(x = Gene, y = AUC)) +
    geom_bar(stat = "identity", fill=paint) +
    coord_flip() +
    ylim(0, 1) +  
    theme_minimal(base_size = 12) +
    labs(
      title = paste(virus_name),
      x = "Gene Symbol",
      y = "AUC"
    ) +
    theme(plot.title = element_text(face = "bold", size = 14))
  
  print(p)
  
  ggsave(file.path(opt$out_dir, paste0("auc_by_",virus_name,".png")), plot = p, width = 7, height = 5, dpi = 300)
}

# Correlation among key Abs (train set)
keep_abs <- intersect(ab_tbl$hid, colnames(merge_train))
if (length(keep_abs) >= 2) {
  cor_mat <- cor(merge_train[, ..keep_abs], method="pearson", use="pairwise.complete.obs")
  colnames(cor_mat) <- ab_tbl$gene[match(colnames(cor_mat), ab_tbl$hid)]
  rownames(cor_mat) <- colnames(cor_mat)
  ragg::agg_png(
    file.path(opt$out_dir, "figure4F_abs_correlation.png"),
    width = 1200, height = 1200, res = 300
  )
  corrplot(cor_mat, method="color", type="upper", tl.cex=1.0, tl.col="black", tl.srt=45, addCoef.col="black")
  dev.off()
}

# ---- main function (drop-in) ----
run_lasso_feature_analysis <- function(
    virus_name,              # "HSV-1" or long name
    antibody_id,             # e.g., "h6808"
    gene_id,                 # e.g., "PHLDA1" (for titles)
    train_vi_anno,           # annotation table with at least UniProt_acc, taxon_species, product
    train_vi, train_hu,      # MGBB/LLF binary matrices with Sample_ID and features
    test_vi, test_hu,        # ABC binary matrices with Sample_ID and features
    out_dir,                 # output directory
    color_map
) {
  message("\nRunning detailed feature LASSO for ", virus_name, " / ", antibody_id, " / ", gene_id)
  
  # standardize species labels in the annotation
  train_vi_anno <- train_vi_anno %>% mutate(taxon_species = species_label_map(taxon_species))
  
  # allow user to pass long or short virus name
  virus_short <- species_label_map(virus_name)
  
  # 1) collect UniProt IDs for this virus present in BOTH train and test
  virus_pro <- train_vi_anno %>%
    filter(taxon_species == virus_short) %>%
    pull(UniProt_acc) %>%
    intersect(colnames(train_vi)) %>%
    intersect(colnames(test_vi))
  
  if (!length(virus_pro)) {
    message("No viral proteins found for ", virus_short, " that are shared in train and test.")
    return(NULL)
  }
  
  # 2) build train & test joined frames
  train_all <- train_hu %>% right_join(train_vi, by = "Sample_ID")
  test_all  <- test_hu  %>% right_join(test_vi,  by = "Sample_ID")
  
  needed_cols <- c(virus_pro, antibody_id)
  if (!all(needed_cols %in% colnames(train_all))) {
    message("Missing columns in training data for ", virus_short, " / ", antibody_id)
    return(NULL)
  }
  if (!all(needed_cols %in% colnames(test_all))) {
    message("Missing columns in test data for ", virus_short, " / ", antibody_id)
    return(NULL)
  }
  
  # keep only required columns; binarize defensively
  df_train <- train_all %>% select(all_of(needed_cols)) %>% binarize01_df()
  df_test  <- test_all  %>% select(all_of(needed_cols))  %>% binarize01_df()
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
  cv_fit <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 1, nfolds = 5, type.measure = "auc")
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
  
  # 7) coefficient plot
  p_coef <- ggplot(keys3, aes(x = reorder(variables, -abs(beta)), y = beta, fill = product)) +
    geom_col(color = "black") +
    geom_text(aes(label = sig_marker,
                  y = ifelse(beta > 0, beta + 0.05, beta - 0.05)),
              size = 3) +
    labs(title = paste0("Prediction model of ", virus_short, " peptides for ", gene_id),
         x = paste(virus_short, "peptides"), y = "LASSO weight", fill = "Product") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 8),
          legend.title = element_text(size=8),
          legend.text = element_text(size=7),
          legend.key.size = unit(0.5, "cm"),
          plot.title = element_text(face = "bold", size = 12)) +
    scale_fill_discrete()
  
  ragg::agg_png(file.path(out_dir, paste0("figure4_", virus_short, "_", antibody_id, "_weights.png")),
                width = 2300, height = 1500, res = 300)
  print(p_coef); dev.off()
  
  # 8) ABC per-variable 2×2 metrics (sensitivity/specificity/AUC_1threshold)
  #    Evaluate *each* selected viral peptide individually on ABC
  results_abc <- lapply(keys3$variables, function(v) {
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
               AUC_1threshold=auc1, P.value=pv, stringsAsFactors = FALSE)
  })
  results_abc <- data.table::rbindlist(results_abc, fill = TRUE, use.names = TRUE)
  results_abc <- results_abc[order(-AUC_1threshold), ]
  results_abc1 <- left_join(results_abc, vi_anno_for_plot[, c("UniProt_acc","product")],
                            by = c("Variable"="UniProt_acc"))
  
  # save the table
  data.table::fwrite(results_abc1,
                     file.path(out_dir, paste0("figure4_", virus_short, "_", antibody_id, "_ABC_metrics.tsv")),
                     sep = "\t")
  
  # 9) ABC scatter (Sensitivity vs Specificity, colored by AUC_1threshold)
  plot_data <- results_abc %>% filter(!is.na(Sensitivity), !is.na(Specificity))
  p_roc <- ggplot(plot_data, aes(x = Sensitivity, y = Specificity,
                                 color = AUC_1threshold, label = Variable)) +
    geom_point(size = 2.8) +
    ggrepel::geom_text_repel(size = 2.6, show.legend = FALSE, max.overlaps = 40) +
    scale_color_viridis_c(option = "plasma", limits = c(0.5, 1), name = expression(AUC)) +
    theme_minimal() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = paste0(virus_short, " peptides for ", gene_id, " in ABC cohort"),
         x = "Sensitivity", y = "Specificity") +
    theme(legend.text = element_text(size = 8),
          plot.title  = element_text(hjust = 0.5, size = 12),
          axis.title  = element_text(size = 10),
          axis.text   = element_text(size = 8))
  
  ragg::agg_png(file.path(out_dir, paste0("figure4_", virus_short, "_", antibody_id, "_ROC.png")),
                width = 1600, height = 1600, res = 300)
  print(p_roc); dev.off()
  
  invisible(list(coef = keys2, abc_metrics = results_abc1))
}

# Run for the main 5 panels
pairs_to_run <- list(
  list(virus = "Human alphaherpesvirus 1", ab = "h6808",  Gene = "PHLDA1"),
  list(virus = "Human betaherpesvirus 5",   ab = "h6329", Gene = "ZNF550"),
  list(virus = "Human betaherpesvirus 5",   ab = "h3143",  Gene = "IQCB1"),
  list(virus = "Human betaherpesvirus 5",   ab = "h8069",  Gene = "DNAJC12"),
  list(virus = "Human gammaherpesvirus 4",   ab = "h8250",  Gene = "P3H4")
)

for (pair in pairs_to_run) {
  run_lasso_feature_analysis(pair$virus, pair$ab, pair$Gene,
                             train_vi_anno = train_vi_anno,
                             train_vi = train_vi, train_hu = train_hu,
                             test_vi = test_vi,   test_hu  = test_hu,
                             out_dir = opt$out_dir,
                             color_map = color_map)
}
