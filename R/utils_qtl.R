suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr)
  library(pROC); library(glmnet); library(stringr)
})

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
