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

species_levels <- c("HSV-1","HSV-2","VZV","CMV","HHV-6A","HHV-6B","HHV-7","EBV","HHV-8")

cap_vals <- function(x, cap) pmin(x, cap)

# collapse (antibody x product) to min P for Manhattan
collapse_minP <- function(dt) {
  dt %>%
    mutate(product1 = paste0(gene_symbol, "_", product)) %>%
    group_by(product1) %>%
    slice_min(order_by = P, with_ties = FALSE) %>%
    ungroup()
}


# Disease categories
dx_categories <- list(
  
  Auto_organ = c(
    "Alopecia_areata", "Autoimmune_hemolytic_anemia",
    "Autoimmune_hepatitis", "Celiac_disease", "Graves_disease", "Guillain_Barre_syndrome",
    "Hashimoto_thyroiditis", "Inflammatory_bowel_disease", "Multiple_sclerosis",
    "Myasthenia_gravis", "Neuromyelitis_optica", "Pemphigus", "Pernicious_anemia",
    "Psoriasis", "Type1_Diabetes", "Vitiligo"
  ),
  
  Auto_system = c(
    "Ankylosing_spondylitis", "Behcet_disease", "Dermatomyositis",
    "Giant_cell_arteritis", "Polymyalgia_rheumatica", "Rheumatoid_arthritis", "Scleroderma",
    "Sjogren_syndrome", "Systemic_lupus_erythematosus"
  ),
  
  Oncology = c(
    "Bladder_cancer", "Breast_cancer", "Cancer", "Colon_Rectum_cancer", "Leukemia",
    "Liver_cancer", "Lung_cancer", "Melanoma", "Non_Hodgkin_Lymphoma",
    "Ovary_cancer", "Pancrease_cancer", "Prostate_cancer", "Stomach_cancer",
    "Uterine_corpus_cancer"
  ),
  
  Cardio = c(
    "Aortic_aneurysm_dissection", "Arterial_embolism_thrombosis", "Cardiomyopathy",
    "Chronic_pulmonary_heart_disease", "Congestive_heart_failure",
    "Coronary_artery_disease", "CVD_excl_stroke", "Hypertensive_heart_disease",
    "Myocardial_infarction", "Peripheral_vascular_disease", "Stroke"
  ),
  
  Neuro = c(
    "Alzheimer", "CognitiveDeficit", "FrontotemporalDementia",
    "LewyBody", "NonSpec_Dementia", "Parkinson", "VascularDementia"
  ),
  
  Respiratory = c(
    "Allergic_rhinitis", "Bronchiectasis", "Chronic_Bronchitis",
    "Coin_lesion_lung_disease", "COPD", "Emphysema", "Interstitial_lung_disease"
  ),
  
  Hepatic = c(
    "Alcohol_liver_disease", "Chronic_liver_disease", "Chronic_viral_hepatitis",
    "Cirrhosis", "NAFLD", "Other_chronic_heptitis"
  ),
  
  Metabolic_Renal = c(
    "Chronic_kidney_disease", "Hypoparathyroidism", "Type2_Diabetes"
  )
)
