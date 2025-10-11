#!/usr/bin/env Rscript
# Figure 4: LASSO prediction analysis
# - Trains per-virus LASSO (glmnet) to predict Human Ab positivity from viral peptides
# - External test ROC (ABC), internal ROC (MGBB)
# - Sensitivity/Specificity (one-threshold) summaries
# - β weight barplots with significance stars
# - Correlation heatmap across key Human Abs
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(glmnet)
  library(pROC)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(corrplot)
  library(stringr)
})

# --------------------------
# CLI
# --------------------------
parser <- ArgumentParser(description = "Figure 4 LASSO prediction pipeline")
parser$add_argument("--train_vi_bin", required=TRUE,
                    help="Training viral peptide binary matrix (rows: Subject_Id, cols: UniProt_acc)")
parser$add_argument("--train_hu_bin", required=TRUE,
                    help="Training human full-length Ab binary matrix (rows: Subject_Id, cols: hXXXX)")
parser$add_argument("--test_vi_bin", required=TRUE,
                    help="Test viral peptide binary matrix (rows: Subject_Id, cols: UniProt_acc)")
parser$add_argument("--test_hu_bin", required=TRUE,
                    help="Test human full-length Ab binary matrix (rows: Subject_Id, cols: hXXXX)")
parser$add_argument("--vi_annotation", required=TRUE,
                    help="VirSIGHT Promax annotation CSV (must include UniProt_acc, taxon_species, product)")
parser$add_argument("--hu_annotation", required=TRUE,
                    help="HuSIGHT FullLength annotation CSV (must include gene_symbol, var_id [hXXXX])")
parser$add_argument("--rep_qtl", required=TRUE,
                    help="Significant QTL file (rep_hsv.sig) with columns: taxon_species, antibody, P, etc.")
parser$add_argument("--outdir", required=TRUE, help="Output directory (tables/plots)")
parser$add_argument("--seed", type="integer", default=123,
                    help="Random seed for CV (default: 123)")
parser$add_argument("--min_auc_plot", type="double", default=0.80,
                    help="Only plot ROC per-Ab where max AUC across viruses >= this (default 0.80)")
parser$add_argument("--abs_of_interest", nargs="+", default=c("h6808","h14835","h3143","h8069","h8250"),
                    help="Human Abs to emphasize in weight and correlation plots (default: PHLDA1, ZNF550, IQCB1, DNAJC12, P3H4)")
args <- parser$parse_args()

dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(args$seed)

# --------------------------
# Read data
# --------------------------
message("Reading inputs...")
vi_anno <- fread(args$vi_annotation)       # VirSIGHT promax (peptide rows)
hu_anno <- fread(args$hu_annotation)       # HuSIGHT FullLength (protein rows)
qtl_sig <- fread(args$rep_qtl)

# Make sure these columns exist
stopifnot(all(c("UniProt_acc","taxon_species") %in% colnames(vi_anno)))
if (!"product" %in% colnames(vi_anno)) {
  vi_anno[, product := NA_character_]
}
stopifnot(all(c("var_id") %in% colnames(hu_anno)))
if (!"gene_symbol" %in% colnames(hu_anno)) {
  hu_anno[, gene_symbol := var_id]
}

train_vi <- fread(args$train_vi_bin)
train_hu <- fread(args$train_hu_bin)
test_vi  <- fread(args$test_vi_bin)
test_hu  <- fread(args$test_hu_bin)

# Normalize ID column name
idcol <- "Subject_Id"
if (!idcol %in% colnames(train_vi)) {
  stop("train_vi_bin must have a 'Subject_Id' column")
}
if (!idcol %in% colnames(train_hu)) {
  stop("train_hu_bin must have a 'Subject_Id' column")
}
if (!idcol %in% colnames(test_vi)) {
  stop("test_vi_bin must have a 'Subject_Id' column")
}
if (!idcol %in% colnames(test_hu)) {
  stop("test_hu_bin must have a 'Subject_Id' column")
}

# --------------------------
# Focus on HSV family species present in QTLs
# --------------------------
species_hsv <- c("Human alphaherpesvirus 1", "Human alphaherpesvirus 2", "Human alphaherpesvirus 3",
                 "Human betaherpesvirus 5", "Human betaherpesvirus 6A", "Human betaherpesvirus 6B",
                 "Human betaherpesvirus 7", "Human gammaherpesvirus 4", "Human gammaherpesvirus 8")

vi_anno <- vi_anno %>% filter(taxon_species %in% species_hsv)
qtl_sig <- qtl_sig %>% filter(taxon_species %in% species_hsv)

# Ensure Hu var_ids exist in matrices
all_abs <- intersect(colnames(train_hu), colnames(test_hu))
if (length(all_abs) == 0) stop("No overlapping Human Ab columns between training and test matrices.")
hu_keep <- union(args$abs_of_interest, unique(qtl_sig$antibody))
hu_keep <- intersect(hu_keep, all_abs)

train_hu1 <- train_hu %>% select(all_of(c(idcol, hu_keep)))
test_hu1  <- test_hu  %>% select(all_of(c(idcol, hu_keep)))

# Ensure VI proteins are in matrices
vi_keep <- intersect(colnames(train_vi), colnames(test_vi))
vi_keep <- intersect(vi_keep, vi_anno$UniProt_acc)

train_vi1 <- train_vi %>% select(all_of(c(idcol, vi_keep)))
test_vi1  <- test_vi  %>% select(all_of(c(idcol, vi_keep)))

# --------------------------
# Helper: fit LASSO per virus → Ab
# --------------------------
species_list <- sort(unique(vi_anno$taxon_species))
results_nested <- list()

for (virus_name in species_list) {
  message(">> Virus: ", virus_name)

  virus_pro <- vi_anno %>%
    filter(taxon_species == virus_name) %>%
    pull(UniProt_acc) %>%
    intersect(vi_keep)

  if (length(virus_pro) == 0) {
    message("   (no peptides present in matrices; skip)")
    results_nested[[virus_name]] <- list()
    next
  }

  # align training
  df_train <- train_hu1 %>%
    right_join(train_vi1 %>% select(all_of(c(idcol, virus_pro))), by = idcol)

  # align test
  df_test  <- test_hu1 %>%
    right_join(test_vi1 %>% select(all_of(c(idcol, virus_pro))),  by = idcol)

  ab_results <- list()

  for (ab in hu_keep) {
    if (!ab %in% colnames(df_train) || !ab %in% colnames(df_test)) next
    if (length(unique(df_train[[ab]])) < 2 || length(unique(df_test[[ab]])) < 2) next

    # Build model matrices
    train_x <- model.matrix(~ ., df_train[, virus_pro, with = FALSE])
    train_y <- df_train[[ab]]
    test_x  <- model.matrix(~ ., df_test[,  virus_pro, with = FALSE])
    test_y  <- df_test[[ab]]

    set.seed(args$seed)
    cv_fit <- cv.glmnet(train_x, train_y, family = "binomial", alpha = 1)

    # External test ROC
    pred_prob_test <- as.vector(predict(cv_fit, newx = test_x, s = cv_fit$lambda.min, type = "response"))
    roc_test <- tryCatch(roc(test_y, pred_prob_test, quiet = TRUE), error = function(e) NULL)

    # Internal (train) ROC
    pred_prob_train <- as.vector(predict(cv_fit, newx = train_x, s = cv_fit$lambda.min, type = "response"))
    roc_train <- tryCatch(roc(train_y, pred_prob_train, quiet = TRUE), error = function(e) NULL)

    # Save
    ab_results[[ab]] <- list(
      cv = cv_fit,
      virus = virus_name,
      AUC_test  = if (!is.null(roc_test))  as.numeric(auc(roc_test))  else NA_real_,
      AUC_train = if (!is.null(roc_train)) as.numeric(auc(roc_train)) else NA_real_,
      ROC_test  = roc_test,
      ROC_train = roc_train,
      features  = virus_pro
    )
  }

  results_nested[[virus_name]] <- ab_results
}

# --------------------------
# Collate AUC tables
# --------------------------
auc_summary <- do.call(rbind, lapply(names(results_nested), function(v) {
  res <- results_nested[[v]]
  if (length(res) == 0) return(NULL)
  data.frame(
    Virus = v,
    Antibody = names(res),
    AUC_test  = sapply(res, function(x) x$AUC_test),
    AUC_train = sapply(res, function(x) x$AUC_train)
  )
}))
if (is.null(auc_summary)) stop("No models were successfully fit.")

auc_ann <- auc_summary %>%
  left_join(hu_anno %>% select(var_id, gene_symbol), by = c("Antibody" = "var_id"))

fwrite(auc_ann, file.path(args$outdir, "auc_hsv_summary.tsv"), sep = "\t")

# --------------------------
# Export ROC coordinates (External test)
# --------------------------
roc_by_ab <- list()
for (virus in names(results_nested)) {
  for (ab in names(results_nested[[virus]])) {
    rr <- results_nested[[virus]][[ab]]
    if (is.null(rr$ROC_test)) next
    cd <- coords(rr$ROC_test, x = "all", ret = c("specificity", "sensitivity"))
    roc_df <- as.data.frame(cd)
    roc_df$FPR <- 1 - roc_df$specificity
    roc_df$TPR <- roc_df$sensitivity
    roc_df$Antibody <- ab
    roc_df$Virus <- virus
    roc_df$AUC <- rr$AUC_test
    roc_by_ab[[paste0(virus,"__",ab)]] <- roc_df
  }
}
roc_all <- bind_rows(roc_by_ab)
roc_all <- roc_all %>% left_join(hu_anno %>% select(var_id, gene_symbol), by = c("Antibody" = "var_id"))
fwrite(roc_all, file.path(args$outdir, "roc_hsv_summary_test.tsv"), sep = "\t")

# --------------------------
# Plot ROC per Ab across viruses (External)
# --------------------------
roc_all$Virus <- factor(roc_all$Virus,
  levels=c("Human alphaherpesvirus 1","Human alphaherpesvirus 2","Human alphaherpesvirus 3",
           "Human betaherpesvirus 5","Human betaherpesvirus 6A","Human betaherpesvirus 6B",
           "Human betaherpesvirus 7","Human gammaherpesvirus 4","Human gammaherpesvirus 8"))

best_auc <- auc_ann %>% group_by(Antibody) %>% summarise(maxAUC = max(AUC_test, na.rm=TRUE), .groups="drop")
to_plot_abs <- best_auc %>% filter(maxAUC >= args$min_auc_plot) %>% pull(Antibody)
if (length(to_plot_abs) == 0) to_plot_abs <- unique(auc_ann$Antibody)

for (ab in to_plot_abs) {
  sub <- roc_all %>% filter(Antibody == ab)
  if (nrow(sub) == 0) next
  gsym <- unique(sub$gene_symbol)
  p <- ggplot(sub, aes(x = FPR, y = TPR, color = Virus)) +
    geom_line(linewidth = 1) +
    geom_abline(linetype = "dashed", color = "gray20") +
    theme_minimal(base_size = 12) +
    labs(title = paste0("ROC (External) for ", paste(gsym, collapse = ","), " — max AUC: ", sprintf("%.3f", max(sub$AUC, na.rm = TRUE))),
         x = "False Positive Rate", y = "True Positive Rate") +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "right")
  ggsave(file.path(args$outdir, paste0("roc_external_by_", ab, ".png")), p, width = 5.5, height = 5, dpi = 300)
}

# --------------------------
# Correlation (train) among key Abs
# --------------------------
abs_for_corr <- intersect(args$abs_of_interest, colnames(train_hu1))
if (length(abs_for_corr) >= 2) {
  cm <- suppressWarnings(cor(train_hu1[, abs_for_corr, with=FALSE], method = "pearson"))
  colnames(cm) <- hu_anno$gene_symbol[match(colnames(cm), hu_anno$var_id)]
  rownames(cm) <- colnames(cm)
  png(file.path(args$outdir, "corr_abs_train.png"), width = 1100, height = 900, res = 150)
  corrplot(cm, method = "color", type = "upper", tl.cex = 0.8, tl.col = "black",
           tl.srt = 45, addCoef.col = "black")
  dev.off()
}

# ==============================================================================
# β weight barplots + one-threshold sensitivity/specificity summaries
#   for highlighted Abs (default: PHLDA1, ZNF550, IQCB1, DNAJC12, P3H4)
#   If a virus is not specified, pick the virus with best external AUC for that Ab.
# ==============================================================================

# Helper to clean product labels (like in your script)
clean_product <- function(x) {
  x <- str_replace_all(x, fixed("Alkaline nuclease (EC 3.1.-.-)"), "Alkaline nuclease")
  x <- str_replace_all(x, fixed("Capsid assembly protein UL37"), "Capsid assembly protein")
  x <- str_replace_all(x, fixed("Deneddylase (EC 3.4.19.12)"), "Deneddylase")
  x <- str_replace_all(x, fixed("E3 ubiquitin-protein ligase ICP0 (EC 6.3.2.-)"), "E3 ubiquitin-protein ligase")
  x <- str_replace_all(x, fixed("Envelope protein UL45 (18 kDa protein)"), "Envelope protein")
  x <- str_replace_all(x, fixed("Envelope glycoprotein E"), "Envelope glycoprotein")
  x <- str_replace_all(x, fixed("Envelope glycoprotein I"), "Envelope glycoprotein")
  x <- str_replace_all(x, fixed("Large tegument protein deneddylase (EC 3.4.19.12) (EC 3.4.22.-)"), "Large tegument protein deneddylase")
  x <- str_replace_all(x, fixed("Large tegument protein deneddylase;Large tegument protein"), "Large tegument protein deneddylase")
  x <- str_replace_all(x, fixed("Neurovirulence factor ICP34.5 (Infected cell protein 34.5) (protein gamma(1)34.5)"), "Neurovirulence protein")
  x <- str_replace_all(x, fixed("Neurovirulence protein ICP34.5"), "Neurovirulence protein")
  x <- str_replace_all(x, fixed("Thymidine kinase (EC 2.7.1.21)"), "Thymidine kinase")
  x <- str_replace_all(x, fixed("US11 *1;US11 *1"), "US11")
  x <- str_replace_all(x, fixed("Uracil-DNA glycosylase (UDG) (EC 3.2.2.27) (UNG)"), "Uracil-DNA glycosylase")
  x <- str_replace_all(x, fixed("Major DNA-binding protein;UL29"), "Major DNA-binding protein")
  x
}

# Choose “best virus per Ab” from external AUC
best_map <- auc_ann %>%
  group_by(Antibody) %>%
  slice_max(order_by = AUC_test, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Antibody, Virus)

# For manuscript defaults: force mapping if present in results
force_map <- tribble(
  ~Antibody, ~Virus,
  "h6808",  "Human alphaherpesvirus 1", # PHLDA1 ~ HSV-1
  "h14835", "Human betaherpesvirus 5",  # ZNF550 ~ CMV
  "h3143",  "Human betaherpesvirus 5",  # IQCB1  ~ CMV
  "h8069",  "Human betaherpesvirus 5",  # DNAJC12 ~ CMV
  "h8250",  "Human gammaherpesvirus 4"  # P3H4   ~ EBV
)
# use forced mapping if that virus exists in results; otherwise fallback to best_map
map_use <- left_join(force_map, best_map, by="Antibody", suffix=c("_force","_best")) %>%
  transmute(Antibody,
            Virus = ifelse(!is.na(Virus_force), Virus_force, Virus_best)) %>%
  distinct()

# helper: one-threshold metrics (2x2 table) in a given cohort df
one_threshold_metrics <- function(df, outcome, marker) {
  out <- data.frame(Variable = marker,
                    Positive_number = NA_real_,
                    Shared_number   = NA_real_,
                    Sensitivity     = NA_real_,
                    Specificity     = NA_real_,
                    AUC_1threshold  = NA_real_,
                    P.value         = NA_real_)
  if (!all(c(outcome, marker) %in% colnames(df))) return(out)

  x <- as.factor(df[[marker]])
  y <- as.factor(df[[outcome]])
  tb <- table(y, x)
  if (!all(dim(tb) == c(2,2))) return(out)

  test1 <- suppressWarnings(chisq.test(tb))
  out$P.value <- as.numeric(test1$p.value)

  TN <- as.numeric(tb[1,1]); FP <- as.numeric(tb[1,2])
  FN <- as.numeric(tb[2,1]); TP <- as.numeric(tb[2,2])

  out$Sensitivity <- ifelse((TP+FN) > 0, TP/(TP+FN), NA_real_)
  out$Specificity <- ifelse((TN+FP) > 0, TN/(TN+FP), NA_real_)
  out$AUC_1threshold <- mean(c(out$Sensitivity, out$Specificity), na.rm=TRUE)
  out$Shared_number <- TP
  out$Positive_number <- TP + FP
  out
}

# produce β barplot + one-threshold scatter for a (virus, antibody)
make_weight_and_threshold_plots <- function(virus_name, ab) {
  rr <- results_nested[[virus_name]][[ab]]
  if (is.null(rr)) return(invisible(NULL))

  cv_fit <- rr$cv
  feats  <- rr$features

  # Extract non-zero betas at lambda.min
  co <- as.matrix(coef(cv_fit, s = "lambda.min"))
  co <- data.frame(beta = as.numeric(co), variables = rownames(co), row.names = NULL)
  co <- co %>% filter(beta != 0, variables != "(Intercept)")

  # Join product annotations
  keys1 <- co %>%
    inner_join(vi_anno %>% select(UniProt_acc, taxon_species, product),
               by = c("variables" = "UniProt_acc"))

  keys1$product <- factor(clean_product(keys1$product))

  # OLS p-values on selected predictors (approx for display)
  # Build train design with selected vars only
  virus_pro <- feats
  df_train <- train_hu1 %>%
    right_join(train_vi1 %>% select(all_of(c(idcol, virus_pro))), by = idcol)
  selected_vars <- keys1$variables
  if (length(selected_vars) > 0) {
    Xsel <- model.matrix(reformulate(termlabels = selected_vars), data = df_train)
    ysel <- df_train[[ab]]
    # Use linear regression for p-value display (as in your script)
    ols_model <- lm(ysel ~ Xsel)
    sm <- summary(ols_model)
    pv <- data.frame(variables = rownames(sm$coefficients),
                     P_Value  = sm$coefficients[,4], row.names = NULL)
    # row names like "Xselvar1" — pull original names
    pv$variables <- str_remove(pv$variables, "^Xsel")
    keys2 <- keys1 %>%
      left_join(pv, by = "variables") %>%
      mutate(sig_marker = case_when(
        P_Value < 0.001 ~ "***",
        P_Value < 0.01  ~ "**",
        P_Value < 0.05  ~ "*",
        TRUE ~ ""
      ))
  } else {
    keys2 <- keys1 %>% mutate(P_Value = NA_real_, sig_marker = "")
  }

  # β barplot
  title_gene <- hu_anno$gene_symbol[match(ab, hu_anno$var_id)]
  p_beta <- ggplot(keys2, aes(x = reorder(variables, -abs(beta)), y = beta, fill = product)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sig_marker,
                  y = ifelse(beta > 0, beta + 0.05, beta - 0.05)),
              size = 4, color = "black", vjust = ifelse(keys2$beta>0, 0, 1)) +
    labs(title = paste0("LASSO weights for ", title_gene, " (", ab, ") — ", virus_name),
         x = paste0(gsub("^Human ", "", virus_name), " peptides"),
         y = "Weight (β)", fill = "Product") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          plot.title = element_text(hjust = 0.5, size = 14))
  ggsave(file.path(args$outdir, paste0("weights_", ab, "_", gsub("\\s+","_",virus_name), ".png")),
         p_beta, width = 9, height = 5.5, dpi = 300)

  # One-threshold (ABC external)
  # For each selected peptide marker, compute sens/spec vs Ab in TEST set
  df_test <- test_hu1 %>%
    right_join(test_vi1 %>% select(all_of(c(idcol, feats))), by = idcol)

  one_thr_list <- lapply(selected_vars, function(v) one_threshold_metrics(df_test, outcome = ab, marker = v))
  one_thr <- bind_rows(one_thr_list)
  one_thr <- one_thr %>% left_join(vi_anno %>% select(UniProt_acc, product), by = c("Variable" = "UniProt_acc"))
  fwrite(one_thr, file.path(args$outdir, paste0("one_threshold_", ab, "_", gsub("\\s+","_",virus_name), "_ABC.tsv")), sep = "\t")

  plot_data <- one_thr %>% filter(!is.na(Sensitivity), !is.na(Specificity))
  if (nrow(plot_data) > 0) {
    p_sc <- ggplot(plot_data, aes(x = Sensitivity, y = Specificity, color = AUC_1threshold, label = Variable)) +
      geom_point(size = 3) +
      geom_text_repel(size = 3, show.legend = FALSE) +
      scale_color_viridis_c(option = "plasma", limits = c(0.6, 1), name = expression(AUC)) +
      theme_minimal() +
      coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
      labs(title = paste0(gsub("^Human ", "", virus_name), " peptides for ", title_gene, " (ABC cohort)"),
           x = "Sensitivity", y = "Specificity")
    ggsave(file.path(args$outdir, paste0("one_threshold_scatter_", ab, "_", gsub("\\s+","_",virus_name), "_ABC.png")),
           p_sc, width = 7, height = 5.5, dpi = 300)
  }

  # One-threshold (MGBB internal)
  df_train2 <- train_hu1 %>%
    right_join(train_vi1 %>% select(all_of(c(idcol, feats))), by = idcol)
  one_thr_list2 <- lapply(selected_vars, function(v) one_threshold_metrics(df_train2, outcome = ab, marker = v))
  one_thr2 <- bind_rows(one_thr_list2)
  one_thr2 <- one_thr2 %>% left_join(vi_anno %>% select(UniProt_acc, product), by = c("Variable" = "UniProt_acc"))
  fwrite(one_thr2, file.path(args$outdir, paste0("one_threshold_", ab, "_", gsub("\\s+","_",virus_name), "_MGBB.tsv")), sep = "\t")

  plot_data2 <- one_thr2 %>% filter(!is.na(Sensitivity), !is.na(Specificity))
  if (nrow(plot_data2) > 0) {
    p_sc2 <- ggplot(plot_data2, aes(x = Sensitivity, y = Specificity, color = AUC_1threshold, label = Variable)) +
      geom_point(size = 3) +
      geom_text_repel(size = 3, show.legend = FALSE) +
      scale_color_gradientn(colours = c("blue","green","yellow","red"), name = expression(AUC)) +
      theme_minimal() +
      coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
      labs(title = paste0(gsub("^Human ", "", virus_name), " peptides for ", title_gene, " (MGBB cohort)"),
           x = "Sensitivity", y = "Specificity")
    ggsave(file.path(args$outdir, paste0("one_threshold_scatter_", ab, "_", gsub("\\s+","_",virus_name), "_MGBB.png")),
           p_sc2, width = 7, height = 5.5, dpi = 300)
  }
}

# Drive β/threshold plots for each Ab using mapping
for (ab in intersect(args$abs_of_interest, names(best_map$Antibody %>% {auc_ann$Antibody}))) {
  vv <- map_use$Virus[match(ab, map_use$Antibody)]
  if (is.na(vv)) {
    # fallback to best virus by AUC
    vv <- best_map$Virus[match(ab, best_map$Antibody)]
  }
  if (!is.na(vv) && vv %in% names(results_nested) && ab %in% names(results_nested[[vv]])) {
    make_weight_and_threshold_plots(vv, ab)
  }
}

# --------------------------
# Save session info for reproducibility
# --------------------------
writeLines(c(capture.output(sessionInfo())),
           con = file.path(args$outdir, "sessionInfo_figure4.txt"))

message("Done. Outputs written to: ", args$outdir)

