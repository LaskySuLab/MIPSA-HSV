#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(dplyr); library(glmnet)
  library(pROC); library(ggplot2); library(ggrepel); library(scales); library(corrplot); library(purrr); library(stringr)
})

option_list <- list(
  make_option("--train_vi", type="character", help="MGBB virus_proteins_binary.txt"),
  make_option("--train_hu", type="character", help="MGBB human_fl_binary.txt"),
  make_option("--train_vi_anno", type="character", help="MGBB VirSIGHT FoB CSV"),
  make_option("--train_hu_anno", type="character", help="MGBB HuSIGHT FL CSV"),
  make_option("--test_vi", type="character", help="ABC virus binary matrix"),
  make_option("--test_hu", type="character", help="ABC human binary matrix"),
  make_option("--ab_list", type="character",
              help="Comma sep gene:hid pairs, e.g. 'PHLDA1:6808,ZNF550:14835,IQCB1:3143,DNAJC12:8069,P3H4:8250'"),
  make_option("--species_map", type="character",
              help="Semicolon sep 'Human alphaherpesvirus 1=HSV-1;Human betaherpesvirus 5=CMV;Human gammaherpesvirus 4=EBV'"),
  make_option("--out_dir", type="character", default="results/Figure4", help="Output dir")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$out_dir, showWarnings=FALSE, recursive=TRUE)

train_vi <- fread(opt$train_vi); train_hu <- fread(opt$train_hu)
test_vi  <- fread(opt$test_vi);  test_hu  <- fread(opt$test_hu)

train_hu_anno <- fread(opt$train_hu_anno)
train_vi_anno <- fread(opt$train_vi_anno)

if (!"var_id" %in% names(train_hu_anno)) train_hu_anno[, var_id := paste0("h", seq_len(.N))]

# Map long names to short tags for facet/legend
species_map <- str_split(opt$species_map, ";")[[1]] |>
  setNames(nm = sub("=.*$", "", .)) |>
  lapply(\(x) sub("^.*=", "", x)) |> unlist()

# Focus antibodies
ab_pairs <- str_split(opt$ab_list, ",")[[1]]
ab_tbl   <- tibble::tibble(
  gene = sub(":.*$", "", ab_pairs),
  hid  = paste0("h", sub("^.*:", "", ab_pairs))
)

# Shared subjects join for safe modeling
merge_train <- train_hu %>%
  right_join(train_vi, by="Subject_Id") %>% as.data.table()
merge_test  <- test_hu %>% right_join(test_vi, by="Subject_Id") %>% as.data.table()

# virus species sets per short label
species_sets <- lapply(names(species_map), function(sp) {
  acc <- train_vi_anno %>% filter(taxon_species == sp) %>% pull(UniProt_acc) %>% intersect(colnames(train_vi))
  tibble::tibble(species = sp, short = species_map[[sp]], acc = list(acc))
}) %>% bind_rows()

# A helper to fit lasso and compute ROC on test
fit_lasso_eval <- function(species_acc, ab_id, df_train, df_test, out_prefix) {
  acc <- species_acc
  if (length(acc) == 0) return(NULL)
  if (!(ab_id %in% colnames(df_train)) || !(ab_id %in% colnames(df_test))) return(NULL)
  if (length(unique(df_train[[ab_id]])) < 2 || length(unique(df_test[[ab_id]])) < 2) return(NULL)

  train_x <- model.matrix(~ ., df_train[, ..acc])
  train_y <- df_train[[ab_id]]
  test_x  <- model.matrix(~ ., df_test[,  ..acc])
  test_y  <- df_test[[ab_id]]

  set.seed(123)
  cv_fit <- cv.glmnet(train_x, train_y, family="binomial", alpha=1)
  pred_prob <- as.numeric(predict(cv_fit, newx=test_x, s=cv_fit$lambda.min, type="response"))

  roc_obj <- tryCatch(roc(test_y, pred_prob, quiet = TRUE), error=function(e) NULL)
  if (is.null(roc_obj)) return(NULL)

  auc_val <- as.numeric(auc(roc_obj))
  roc_df  <- coords(roc_obj, x="all", ret=c("specificity","sensitivity")) %>%
    as.data.frame() %>% mutate(FPR = 1 - specificity, TPR = sensitivity)

  list(roc = roc_df, auc = auc_val, lambda = cv_fit$lambda.min, cv = cv_fit)
}

# Run per species × antibody and collect outputs
roc_all <- list()
auc_tab <- list()

for (i in seq_len(nrow(species_sets))) {
  sp  <- species_sets$species[i]
  ssp <- species_sets$short[i]
  acc <- species_sets$acc[[i]]

  for (j in seq_len(nrow(ab_tbl))) {
    ab  <- ab_tbl$hid[j]
    gene <- ab_tbl$gene[j]

    res <- fit_lasso_eval(acc, ab, merge_train, merge_test, paste0(ssp,"_",ab))
    if (is.null(res)) next

    roc_all[[paste0(ab,"_",ssp)]] <- res$roc %>% mutate(Antibody=ab, Gene=gene, Virus=ssp, AUC=res$auc)
    auc_tab[[paste0(ab,"_",ssp)]] <- data.frame(Antibody=ab, Gene=gene, Virus=ssp, AUC=res$auc)
  }
}

roc_all_df <- bind_rows(roc_all)
auc_tab_df <- bind_rows(auc_tab)

fwrite(auc_tab_df, file.path(opt$out_dir, "auc_summary_ABCtest.csv"))

# ROC plots per antibody with lines per virus
if (nrow(roc_all_df)) {
  for (ab in unique(roc_all_df$Antibody)) {
    sub <- roc_all_df %>% filter(Antibody==ab)
    gene <- unique(sub$Gene)
    p <- ggplot(sub, aes(FPR, TPR, color = Virus)) +
      geom_line(size=1) + geom_abline(linetype="dashed") +
      labs(title = paste0("ROC for ", gene, " (ABC test)  MaxAUC=", sprintf("%.3f", max(sub$AUC))),
           x="False Positive Rate", y="True Positive Rate") +
      theme_minimal(base_size=12)
    ggsave(file.path(opt$out_dir, paste0("roc_", ab, ".png")), p, width=5.5, height=5, dpi=300)
  }
}

# Correlation among key Abs (train set)
keep_abs <- intersect(ab_tbl$hid, colnames(merge_train))
if (length(keep_abs) >= 2) {
  cor_mat <- cor(merge_train[, ..keep_abs], method="pearson", use="pairwise.complete.obs")
  colnames(cor_mat) <- ab_tbl$gene[match(colnames(cor_mat), ab_tbl$hid)]
  rownames(cor_mat) <- colnames(cor_mat)
  png(file.path(opt$out_dir, "figure4F_abs_correlation.png"), width=900, height=900, res=150)
  corrplot(cor_mat, method="color", type="upper", tl.cex=0.9, tl.col="black", tl.srt=45, addCoef.col="black")
  dev.off()
}

# One-threshold metrics per selected variable (sensitivity, specificity, AUC_1threshold)
# For ABC test, compute per-variable 2x2 against Ab
compute_metrics <- function(df, ab, vars) {
  purrr::map_dfr(vars, function(v) {
    if (!v %in% colnames(df)) return(NULL)
    tt <- table(df[[ab]], df[[v]])
    if (!all(dim(tt)==c(2,2))) return(NULL)
    TN <- as.numeric(tt[1,1]); FP <- as.numeric(tt[1,2])
    FN <- as.numeric(tt[2,1]); TP <- as.numeric(tt[2,2])
    sens <- TP/(TP+FN); spec <- TN/(TN+FP)
    data.frame(Variable=v, Sensitivity=sens, Specificity=spec, AUC_1threshold=(sens+spec)/2)
  })
}

# For each ab × species, take non-zero lasso features (by 1se or min) from train fit and evaluate per-variable metrics on ABC
metrics_list <- list()
for (i in seq_len(nrow(species_sets))) {
  ssp <- species_sets$short[i]
  acc <- species_sets$acc[[i]]
  for (j in seq_len(nrow(ab_tbl))) {
    ab <- ab_tbl$hid[j]
    key <- paste0(ab,"_",ssp)
    roc_key <- names(roc_all)[startsWith(names(roc_all), key)]
    if (!length(roc_key)) next
    # refit to retrieve coefs (train)
    tx <- model.matrix(~ ., merge_train[, ..acc])
    ty <- merge_train[[ab]]
    set.seed(123)
    cv_fit <- cv.glmnet(tx, ty, family="binomial", alpha=1)
    coef_min <- as.matrix(coef(cv_fit, s = cv_fit$lambda.min))
    keep <- rownames(coef_min)[as.numeric(coef_min)!=0]
    keep <- setdiff(keep, "(Intercept)")
    if (!length(keep)) next
    met <- compute_metrics(merge_test, ab, intersect(keep, colnames(merge_test)))
    if (!is.null(met) && nrow(met)) {
      met$Antibody <- ab; met$Gene <- ab_tbl$gene[j]; met$Virus <- ssp
      metrics_list[[length(metrics_list)+1]] <- met
    }
  }
}
if (length(metrics_list)) {
  metrics_df <- bind_rows(metrics_list)
  fwrite(metrics_df, file.path(opt$out_dir, "abc_variable_metrics.csv"))
}
