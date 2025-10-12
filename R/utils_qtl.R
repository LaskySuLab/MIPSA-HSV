suppressPackageStartupMessages({
  library(data.table); library(dplyr)
  library(stats); library(purrr)
})

# Build combined disease indicator from *_at_collect and *_flag columns
build_disease_combine <- function(phe_dt) {
  status_vars   <- grep("_Status$", colnames(phe_dt), value = TRUE)
  collect_vars  <- grep("_at_collect$", colnames(phe_dt), value = TRUE)
  disease       <- gsub("_Status$", "", status_vars)
  for (prefix in disease) {
    matched <- c(paste0(prefix,"_at_collect"), paste0(prefix,"_flag"))
    matched <- matched[matched %in% names(phe_dt)]
    if (length(matched)) {
      phe_dt[, (paste0(prefix,"_combine")) :=
               as.integer(rowSums(.SD, na.rm = TRUE) > 0), .SDcols = matched]
    }
  }
  phe_dt
}

recode_covariates <- function(dt) {
  dt <- dt %>%
    mutate(
      Race_Group = ifelse(Race_Group == "White", 0L, 1L),
      Gender     = ifelse(Gender == "M", 0L, 1L)
    )
  # Smoking
  smoking <- c("Current-smoker", "Former-smoker")
  dt <- as.data.table(dt)
  dt[, Smoking_flag_Biobank := ifelse(Smoking_flag_Biobank %in% smoking, "Current/Former-smoker", Smoking_flag_Biobank)]
  dt[Smoking_flag_Biobank == "", Smoking_flag_Biobank := NA]
  dt[, Smoking_RPDR_most_recent := ifelse(Smoking_RPDR_most_recent %in% smoking, "Current/Former-smoker", Smoking_RPDR_most_recent)]
  dt[Smoking_RPDR_most_recent == "", Smoking_RPDR_most_recent := NA]
  dt[, Smoking_flag := Smoking_flag_Biobank]
  dt[is.na(Smoking_flag), Smoking_flag := Smoking_RPDR_most_recent]
  dt[is.na(Smoking_flag), Smoking_flag := "Unknown"]

  # Alcohol
  alcohol.2 <- c("Current-drinker", "Former-drinker")
  noalcohol <- " 1. None, or less than 1 per month "
  dt[, Alcohol_flag_Biobank := ifelse(Alcohol_flag_Biobank %in% alcohol.2, "Current/Former-drinker", Alcohol_flag_Biobank)]
  dt[, Alcohol_flag_Biobank := ifelse(Alcohol_flag_Biobank %in% noalcohol, "No-alcohol", Alcohol_flag_Biobank)]
  dt[Alcohol_flag_Biobank == "", Alcohol_flag_Biobank := NA]
  dt[, Alcohol_RPDR_most_recent := ifelse(Alcohol_RPDR_most_recent %in% alcohol.2, "Current/Former-drinker", Alcohol_RPDR_most_recent)]
  dt[Alcohol_RPDR_most_recent == "", Alcohol_RPDR_most_recent := NA]
  dt[, Alcohol_flag := Alcohol_flag_Biobank]
  dt[is.na(Alcohol_flag), Alcohol_flag := Alcohol_RPDR_most_recent]
  dt[is.na(Alcohol_flag), Alcohol_flag := "Unknown"]

  dt[, Smoking := ifelse(Smoking_flag %in% c("Current/Former-smoker"), 1L, 0L)]
  dt[, Alcohol := ifelse(Alcohol_flag %in% c("Current/Former-drinker"), 1L, 0L)]
  dt[]
}

bh_fdr <- function(p) p.adjust(p, method = "BH")

glm_safe <- function(formula, data) {
  out <- tryCatch(glm(formula, family = "binomial", data = data), error = function(e) NULL)
  if (is.null(out)) return(NULL)
  s <- summary(out)
  if (nrow(s$coefficients) == 8) return(NULL)
  cbind(as.data.frame(s$coefficients[2, 1:4]), row.names = NULL) |>
    setNames(c("Beta", "SE", "Z", "P"))
}
