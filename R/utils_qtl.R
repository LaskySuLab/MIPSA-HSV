suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(pROC)
  library(glmnet)
})

HSV_SPECIES <- c(
  "Human alphaherpesvirus 1","Human alphaherpesvirus 2","Human alphaherpesvirus 3",
  "Human betaherpesvirus 5","Human betaherpesvirus 6A","Human betaherpesvirus 6B",
  "Human betaherpesvirus 7","Human gammaherpesvirus 4","Human gammaherpesvirus 8"
)

species_short <- function(x) {
  out <- x
  maps <- c(
    "Human alphaherpesvirus 1"="HSV-1","Human alphaherpesvirus 2"="HSV-2",
    "Human alphaherpesvirus 3"="VZV","Human betaherpesvirus 5"="CMV",
    "Human betaherpesvirus 6A"="HHV-6A","Human betaherpesvirus 6B"="HHV-6B",
    "Human betaherpesvirus 7"="HHV-7","Human gammaherpesvirus 4"="EBV",
    "Human gammaherpesvirus 8"="HHV-8"
  )
  out[x %in% names(maps)] <- maps[x[x %in% names(maps)]]
  out
}

detect_and_unify_id <- function(dt) {
  nm <- names(dt)
  idcol <- if ("Sample_ID" %in% nm) "Sample_ID" else if ("Subject_Id" %in% nm) "Subject_Id" else NA
  if (is.na(idcol)) stop("ID column not found (expected Sample_ID or Subject_Id).")
  setnames(dt, old = idcol, new = "Sample_ID")
  dt
}

recode_covariates <- function(DT) {
  DT$Race_Group <- ifelse(DT$Race_Group == "White", 0, 1)
  DT$Gender     <- ifelse(DT$Gender     == "M",     0, 1)

  smoking <- c("Current-smoker", "Former-smoker")
  DT[, Smoking_flag_Biobank := ifelse(Smoking_flag_Biobank %in% smoking, "Current/Former-smoker", Smoking_flag_Biobank)]
  DT[Smoking_flag_Biobank=="", Smoking_flag_Biobank:=NA]
  DT[, Smoking_RPDR_most_recent := ifelse(Smoking_RPDR_most_recent %in% smoking, "Current/Former-smoker", Smoking_RPDR_most_recent)]
  DT[Smoking_RPDR_most_recent=="", Smoking_RPDR_most_recent:=NA]

  DT[, Smoking_flag := Smoking_flag_Biobank]
  DT[is.na(Smoking_flag), Smoking_flag := Smoking_RPDR_most_recent]
  DT[is.na(Smoking_flag), Smoking_flag := "Unknown"]

  alcohol_full <- sort(unique(DT$Alcohol_flag_Biobank))
  if (length(alcohol_full)) {
    # Matches your original collapsing
    noalcohol <- " 1. None, or less than 1 per month "
    DT[, Alcohol_flag_Biobank := ifelse(Alcohol_flag_Biobank %in% alcohol_full[alcohol_full != "" & !is.na(alcohol_full)][2:min(9, length(alcohol_full))],
                                        "Current/Former-drinker", Alcohol_flag_Biobank)]
    DT[, Alcohol_flag_Biobank := ifelse(Alcohol_flag_Biobank %in% noalcohol, "No-alcohol", Alcohol_flag_Biobank)]
  }
  DT[Alcohol_flag_Biobank=="", Alcohol_flag_Biobank:=NA]
  DT[, Alcohol_RPDR_most_recent := ifelse(Alcohol_RPDR_most_recent %in% c("Current-drinker", "Former-drinker"),
                                          "Current/Former-drinker", Alcohol_RPDR_most_recent)]
  DT[Alcohol_RPDR_most_recent=="", Alcohol_RPDR_most_recent:=NA]

  DT[, Alcohol_flag := Alcohol_flag_Biobank]
  DT[is.na(Alcohol_flag), Alcohol_flag := Alcohol_RPDR_most_recent]
  DT[is.na(Alcohol_flag), Alcohol_flag := "Unknown"]

  DT[, Smoking := ifelse(Smoking_flag %in% "Current/Former-smoker", 1, 0)]
  DT[, Alcohol := ifelse(Alcohol_flag %in% "Current/Former-drinker", 1, 0)]
  DT
}

#' LASSO prediction of autoantigen from a virus's peptide panel
fit_lasso_predict <- function(X_bin, y_vec, nfolds=5, family="binomial") {
  stopifnot(nrow(X_bin) == length(y_vec))
  cv <- cv.glmnet(as.matrix(X_bin), y_vec, family=family, alpha=1, nfolds=nfolds)
  mod <- glmnet(as.matrix(X_bin), y_vec, family=family, alpha=1, lambda=cv$lambda.min)
  list(model=mod, lambda=cv$lambda.min, cv=cv)
}

#' AUC using pROC
auc_prob <- function(prob, y) {
  if (length(unique(y)) < 2) return(NA_real_)
  as.numeric(pROC::roc(y, prob, quiet=TRUE)$auc)
}

#' Safe glm with binomial family
safe_glm <- function(formula, data) {
  out <- tryCatch(glm(formula, family=binomial(), data=data), error=function(e) NULL)
  if (is.null(out)) return(NULL)
  if (nrow(summary(out)$coefficients) < 2) return(NULL)
  out
}
