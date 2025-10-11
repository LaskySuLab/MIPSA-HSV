suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr)
  library(readr); library(stringr); library(ggplot2)
})

#' Read CSV/TSV by extension
read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv"))   return(data.table::fread(path))
  if (ext %in% c("tsv","txt")) return(data.table::fread(path, sep="\t"))
  stop("Unsupported file extension for: ", path)
}

#' Binarize numeric fold-change/hit scores to 0/1 with threshold and min prevalence
binarize_hits <- function(dt, id_col, thresh=1L, min_prev=0.01) {
  stopifnot(id_col %in% names(dt))
  num_cols <- setdiff(names(dt), id_col)
  x <- copy(dt)
  for (j in num_cols) {
    x[[j]] <- as.integer(x[[j]] >= thresh)
  }
  # drop markers < min_prev
  prev <- colMeans(x[, ..num_cols], na.rm=TRUE)
  keep <- names(prev[prev >= min_prev])
  cbind(x[, ..id_col], x[, ..keep])
}

#' Make covariates (harmonize race, sex, smoking, alcohol) as in manuscript
make_covariates <- function(phe,
                            id_col="Sample_ID",
                            race_col="Race_Group",
                            sex_col="Gender",
                            bmi_col="BMI_most_recent",
                            smoke_cols=c("Smoking_flag_Biobank","Smoking_RPDR_most_recent"),
                            alcohol_cols=c("Alcohol_flag_Biobank","Alcohol_RPDR_most_recent")) {

  df <- as.data.table(phe)
  setnames(df, old=id_col, new="Sample_ID", skip_absent=TRUE)

  # Sex/Gender -> 0/1 (M=0, F=1)
  if (sex_col %in% names(df)) df[, Gender01 := ifelse(get(sex_col) %in% c("F","Female"), 1L, 0L)] else df[, Gender01:=NA_integer_]

  # Race -> White=0, Non-White=1
  if (race_col %in% names(df)) df[, Race01 := ifelse(get(race_col) == "White", 0L, 1L)] else df[, Race01:=NA_integer_]

  # Smoking (Current/Former vs other)
  smoking <- c("Current-smoker","Former-smoker")
  s1 <- if (smoke_cols[1] %in% names(df)) df[[smoke_cols[1]]] else NA
  s2 <- if (smoke_cols[2] %in% names(df)) df[[smoke_cols[2]]] else NA
  s <- ifelse(s1 %in% smoking | s2 %in% smoking, "Current/Former-smoker", "Other/Unknown")
  df[, Smoking01 := as.integer(s == "Current/Former-smoker")]

  # Alcohol (Current/Former vs none/unknown)
  alc_now <- c("Current-drinker","Former-drinker")
  a1 <- if (alcohol_cols[1] %in% names(df)) df[[alcohol_cols[1]]] else NA
  a2 <- if (alcohol_cols[2] %in% names(df)) df[[alcohol_cols[2]]] else NA
  a <- ifelse(a1 %in% alc_now | a2 %in% alc_now, "Current/Former-drinker", "Other/Unknown")
  df[, Alcohol01 := as.integer(a == "Current/Former-drinker")]

  # BMI, Age
  if (bmi_col %in% names(df)) df[, BMI := as.numeric(get(bmi_col))] else df[, BMI:=NA_real_]
  if ("Age_at_collect" %in% names(df)) df[, Age := as.numeric(Age_at_collect)] else df[, Age:=NA_real_]

  keep <- c("Sample_ID","Age","Gender01","BMI","Race01","Smoking01","Alcohol01")
  unique(df[, ..keep])
}

#' Build formula string for GLM
glm_formula <- function(y, x) {
  as.formula(paste0(y, " ~ ", x, " + Age + Gender01 + BMI + Race01 + Smoking01 + Alcohol01"))
}

#' Benjamini–Hochberg FDR
bh <- function(p) p.adjust(p, method="BH")
