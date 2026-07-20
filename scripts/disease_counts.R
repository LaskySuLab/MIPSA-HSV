#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse); library(data.table); library(dplyr);
  library(tidyr); library(ggplot2); library(stringr)
})

opts <- list()
opt_list <- list(
  make_option("--llf_phe", type="character",
              help="LLF phenotype CSV (e.g., llf_1289_phe1.tsv)"),
  make_option("--lec_phe", type="character",
              help="LEC phenotype CSV (e.g., lec_763_phe1.tsv)"),
  make_option("--llf_exclude", type="character", default="Depression,Opioid_use_disorder,Alcoholism,toothache,sciatica,Psychosis,loose_tooth,
              Asthma,Headache,Migraine,Post_Traumatic_Stress_Disorder,Bipolar_Disorder,Female_Infertility,Male_Infertility,
              Schizophrenia,Chronic_viral_hepatitis,Other_chronic_heptitis,Alcohol_liver_disease,Hiatus_Hernia,
              Diverticular,Cataract,Presbyopia,Mouth_ulcer,Lyme_Disease",
              help="Comma-separated disease base names to exclude (without suffix) for LLF"),
  make_option("--out_dir", type="character", default="results1/Figure1",
              help="Output directory [default %default]"),
  make_option("--llf_min_total", type="double", default=1289*0.01,
              help="Min total cases (incident+prevalent) to display for LLF [default %default]")
  )

opt <- parse_args(OptionParser(option_list = opt_list))

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

compute_counts <- function(phe_dt, exclude_vec, n_label) {
  # Identify disease base names present in both *_flag and *_at_collect
  flag_cols    <- grep("_flag$",    names(phe_dt), value = TRUE)
  collect_cols <- grep("_at_collect$", names(phe_dt), value = TRUE)
  disease_names <- intersect(str_remove(flag_cols, "_flag$"),
                             str_remove(collect_cols, "_at_collect$"))
  
  # Compute incident / prevalent counts
  case_counts <- lapply(disease_names, function(disease) {
    flag_col    <- paste0(disease, "_flag")
    collect_col <- paste0(disease, "_at_collect")
    data.frame(
      disease         = disease,
      incident_cases  = sum(phe_dt[[flag_col]]    == 1, na.rm = TRUE),
      prevalent_cases = sum(phe_dt[[collect_col]] == 1, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  
  case_counts$incident_percent  <- case_counts$incident_cases  * 100 / n_label
  case_counts$prevalent_percent <- case_counts$prevalent_cases * 100 / n_label
  
  case_counts %>%
    filter(!(disease %in% exclude_vec)) %>%
    mutate(disease = gsub("_", " ", disease)) %>%
    mutate(total_cases = incident_cases + prevalent_cases) %>%
    arrange(desc(total_cases)) %>%
    mutate(disease = factor(disease, levels = unique(disease)))
}

to_long <- function(case_counts_tbl) {
  case_counts_tbl %>%
    select(disease, incident_cases, prevalent_cases) %>%
    pivot_longer(cols = c("incident_cases", "prevalent_cases"),
                 names_to = "case_type", values_to = "count") %>%
    mutate(case_type = recode(case_type,
                              "incident_cases" = "Incident",
                              "prevalent_cases" = "Prevalent"))
}

plot_counts <- function(df_long, title_txt, subtitle_txt, outfile) {
  p <- ggplot(df_long, aes(x = disease, y = count, fill = case_type)) +
    geom_bar(stat = "identity") +
    theme_minimal(base_size = 12) +
    labs(
      # title = title_txt,
      subtitle = subtitle_txt,
      x = "Disease",
      y = "Number of Cases",
      fill = "Case Type"
    ) +
    scale_fill_manual(values = c("Prevalent" = "#1f78b4", "Incident" = "#F8766D")) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, size =8),
      axis.text   = element_text(size = 11),
      axis.title  = element_text(size = 14, face = "bold"),
      plot.title  = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 10), legend.position = c(0.9, 0.5)
    )
  ggsave(outfile, p, width = 10, height = 6.0, bg='white', dpi = 300)
}

# -------------------------------
# LLF (MGBB-LLF)
# -------------------------------
llf <- fread(opt$llf_phe)
if (!"Subject_Id" %in% names(llf)) llf[, Subject_Id := as.character(Subject_Id)]
n_llf <- nrow(llf)

exclude_llf <- strsplit(opt$llf_exclude, ",")[[1]] |> trimws()
counts_llf  <- compute_counts(llf, exclude_llf, n_label = n_llf)

# Keep > threshold
counts_llf_fil <- counts_llf %>% filter(total_cases > opt$llf_min_total)

# Save tables
fwrite(counts_llf,      file.path(opt$out_dir, "llf_case_counts_all.csv"))
fwrite(counts_llf_fil,  file.path(opt$out_dir, "llf_case_counts_filtered.csv"))

# Plot Figure 1C
llf_long <- to_long(counts_llf_fil)
plot_counts(
  llf_long,
  title_txt   = sprintf("%s (N=%d)", opt$llf_title, n_llf),
  subtitle_txt= NULL,
  outfile     = file.path(opt$out_dir, "figure1C_llf_disease_counts.png")
)

# Optional summaries
capture.output(summary(llf$COLLECTION_DATE),
               file = file.path(opt$out_dir, "llf_collection_date_summary.txt"))
capture.output(summary(llf$FU_endDate),
               file = file.path(opt$out_dir, "llf_fu_enddate_summary.txt"))
capture.output(summary(llf$FU_years),
               file = file.path(opt$out_dir, "llf_fu_years_summary.txt"))

# -------------------------------
# LEC (MGBB-LEC)
# -------------------------------
lec <- fread(opt$lec_phe)
n_lec <- nrow(lec)
counts_lec  <- compute_counts(lec, exclude_llf, n_label = n_lec)

# Keep > threshold
counts_lec_fil <- counts_lec %>% filter(disease %in% counts_llf_fil$disease)

# Save tables
fwrite(counts_lec,      file.path(opt$out_dir, "lec_case_counts_all.csv"))
fwrite(counts_lec_fil,  file.path(opt$out_dir, "lec_case_counts_filtered.csv"))

# Plot Figure 1D
lec_long <- to_long(counts_lec_fil)
plot_counts(
  lec_long,
  title_txt   = sprintf("%s (N=%d)", opt$lec_title, (n_lec)),
  subtitle_txt= NULL,
  outfile     = file.path(opt$out_dir, "figure1D_lec_disease_counts.png")
)

# Optional summaries
capture.output(summary(lec$COLLECTION_DATE),
               file = file.path(opt$out_dir, "lec_collection_date_summary.txt"))
capture.output(summary(lec$FU_endDate),
               file = file.path(opt$out_dir, "lec_fu_enddate_summary.txt"))
capture.output(summary(lec$FU_years),
               file = file.path(opt$out_dir, "lec_fu_years_summary.txt"))

message("Done. Wrote figures + tables to: ", normalizePath(opt$out_dir))

#------------------------------------------------------------------------------#
# Table 1
n_llf
# [1] 1289
n_lec
# [1] 763
table(llf$Gender)
# F   M 
# 907 382 
table(lec$Gender)
# F   M 
# 406 357 
mean(llf$Age_at_collect);sd(llf$Age_at_collect)
# [1] 58.44185
# [1] 14.84713
mean(lec$Age_at_collect);sd(lec$Age_at_collect)
# [1] 55.11407
# [1] 16.90126
mean(llf$BMI_most_recent);sd(llf$BMI_most_recent)
# [1] 30.48352
# [1] 7.749385
mean(lec$BMI_most_recent);sd(lec$BMI_most_recent)
# [1] 30.48013
# [1] 6.534815
table(llf$decease_flag)
# 0    1 
# 1074  215 
table(lec$decease_flag)
# 0   1 
# 579 184

table(llf$Race_Group)
# American Indian or Alaska Native                                     Asian                                     Black 
# 3                                        18                                       109 
# Declined Native Hawaiian or Other Pacific Islander                                     Other 
# 6                                         1                                        76 
# Two or More                           Unknown/Missing                                     White 
# 14                                        13                                      1049 
table(lec$Race_Group)
# American Indian or Alaska Native                                     Asian                                     Black 
# 3                                        11                                        75 
# Declined Native Hawaiian or Other Pacific Islander                                     Other 
# 8                                         1                                       101 
# Two or More                           Unknown/Missing                                     White 
# 5                                        14                                       545 

table(llf$Ethnic_Group)
# DECLINED        HISPANIC    Non Hispanic Unknown/Missing 
# 47              22            1190              30 
table(lec$Ethnic_Group)
# DECLINED        HISPANIC    Non Hispanic Unknown/Missing 
# 28              33             686              16 

table(llf.1$Smoking_flag)
# Current/Former-smoker            Non-smoker 
# 564                   725 
table(lec.1$Smoking_flag)
# Current/Former-smoker            Non-smoker               Unknown 
# 254                   505                     4 
table(llf.1$Alcohol_flag)
# Current/Former-drinker             No-alcohol 
# 974                    315 
table(lec.1$Alcohol_flag)
# Current/Former-drinker             No-alcohol 
# 640                    123  
table(llf.1$cci)
# 1   2   3   4   5   6   7 
# 48  81 143 299 345 233 140 
table(lec.1$cci)
# 0   1   2   3   4   5   6   7 
# 62  51  48  97 146 158 136  65 
mean(llf$FU_years);sd(llf$FU_years)
# [1] 8.596255
# [1] 2.954672
mean(lec$FU_years);sd(lec$FU_years)
# [1] 7.856882
# [1] 2.693484
mean(llf.1$ics_trim_totnum_5y);sd(llf.1$ics_trim_totnum_5y)
# [1] 17.89527
# [1] 24.31274

table(llf.1$rx_antiviral_24m)
# 0    1 
# 777 512
table(lec.1$rx_antiviral_24m)
# 0   1 
# 525 238 
table(llf.1$rx_antineoplastic_24m)
# 0    1 
# 987 302 
table(lec.1$rx_antineoplastic_24m)
# 0   1 
# 609 154 
table(llf.1$rx_immunosuppressant_24m)
# 0    1 
# 1055  234
table(lec.1$rx_immunosuppressant_24m)
# 0   1 
# 680  83 

