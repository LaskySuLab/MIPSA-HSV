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

#' Binarize numeric fold-change/hit scores to 0/1 with threshold
binarize01_df <- function(df) {
  out <- as.data.frame(df)
  for (j in seq_along(out)) {
    v <- out[[j]]
    if (is.numeric(v)) {
      v[is.na(v)] <- 0
      v[v != 0] <- 1
      out[[j]] <- v
    }
  }
  out
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

#' Short virus names
species_label_map <- function(x) {
  x <- gsub('Human alphaherpesvirus 1', 'HSV-1', x)
  x <- gsub('Human alphaherpesvirus 2', 'HSV-2', x)
  x <- gsub('Human alphaherpesvirus 3', 'VZV', x)
  x <- gsub('Human betaherpesvirus 5', 'CMV', x)
  x <- gsub('Human betaherpesvirus 6A', 'HHV-6A', x)
  x <- gsub('Human betaherpesvirus 6B', 'HHV-6B', x)
  x <- gsub('Human betaherpesvirus 7', 'HHV-7', x)
  x <- gsub('Human gammaherpesvirus 4', 'EBV', x)
  x <- gsub('Human gammaherpesvirus 8', 'HHV-8', x)
  x
}

#' HSV species and color matching
species_list <- c("HSV-1","HSV-2","VZV","CMV","HHV-6A","HHV-6B","HHV-7","EBV","HHV-8")
color = c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F", "#00B9E3", "#619CFF", "#DB72FB", "#FF61C3")
color_map <- setNames(
  c("#F8766D","#D39200","#93AA00","#00BA38","#00C19F","#00B9E3","#619CFF","#DB72FB","#FF61C3"),
  c("HSV-1","HSV-2","VZV","CMV","HHV-6A","HHV-6B","HHV-7","EBV","HHV-8")
)


#' optional cleanups of product names
clean_replace <- c(
    "Alkaline nuclease (EC 3.1.-.-)" = "Alkaline nuclease",
    "Deneddylase (EC 3.4.19.12)" = "Deneddylase",
    "E3 ubiquitin-protein ligase ICP0 (EC 6.3.2.-)" = "E3 ubiquitin-protein ligase",
    "Large tegument protein deneddylase (EC 3.4.19.12) (EC 3.4.22.-)" = "Large tegument protein deneddylase",
    "Nuclear egress protein 2" = "Nuclear egress protein",
    "Thymidine kinase (EC 2.7.1.21)" = "Thymidine kinase",
    "US11 *1;US11 *1" = "Protein US11",
    "Uracil-DNA glycosylase (UDG) (EC 3.2.2.27) (UNG)" = "Uracil-DNA glycosylase",
    
    "Envelope glycoprotein UL130" = "Envelope glycoprotein",
    "Membrane RL1 protein2 (Membrane protein RL12)" = "Membrane protein",
    "Membrane glycoprotein US3;Membrane glycoprotein US3" = "Membrane glycoprotein",
    "Membrane glycoprotein US7;Membrane glycoprotein US7" = "Membrane glycoprotein",
    "Membrane protein RL12;Membrane protein RL12" = "Membrane protein",
    "Membrane protein RL13;Membrane protein RL13" = "Membrane protein",
    "Membrane protein RL12" = "Membrane protein",
    "Membrane protein UL20" = "Membrane protein",
    "Tegument protein UL43" = "Tegument protein",
    "UL37" = "Protein UL37",
    " immediate early glycoprotein" = "",
    
    "Envelope glycoprotein gpUL55" = "Envelope glycoprotein",
    "Immediate early protein " = "",
    "Membrane RL1 protein2 (Membrane protein RL12)" = "Membrane protein",
    "Membrane RL1 protein2;Membrane RL1 protein2" = "Membrane protein",
    "Membrane glycoprotein US7;Membrane glycoprotein US7" = "Membrane glycoprotein",
    "ORFS326C;ORFS326C;ORFS326C" = "ORFS326C",
    "Transmembrane protein HWLF3" = "Transmembrane protein",
    "UL74 protein" = "Protein UL74",
    
    "BZLF1 (BZLF1 protein)" = "Transcriptional activator",
    "EBNA3C (EBNA3C nuclear protein)" = "EBV nuclear antigen",
    "Nuclear antigen EBNA-3C;Nuclear antigen EBNA-3C" = "EBV nuclear antigen",
    "Nuclear protein EBNA-2" = "EBV nuclear antigen",
    "Protein BRRF1;Protein BRRF1;Protein BRRF1" = "Transcriptional activator",
    "Uncharacterized protein BLRF3" = "Uncharacterized protein"
)
