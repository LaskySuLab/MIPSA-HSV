#!/usr/bin/env Rscript
# Prevalence and Incidence of Diseases in both LLF (MGBB-LLF) and ABC cohorts
# Produces Figure 1C (LLF) and Figure 1D (ABC) + CSV tables

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

opts <- list()
opt_list <- list(
  make_option("--llf_phe", type="character",
              help="LLF phenotype CSV (e.g., MIPSA_Asthma_1290.csv)"),
  make_option("--abc_phe", type="character",
              help="ABC phenotype CSV (e.g., Ab_pheno.csv)"),
  make_option("--out_dir", type="character", default="results",
              help="Output directory [default %default]"),
  make_option("--llf_exclude", type="character", default="AIDS,Depression,HepatitisB,HepatitisC,Migraine,Opioid_use_disorder,Headache,Asthma,Chronic_viral_hepatitis",
              help="Comma-separated disease base names to exclude (without suffix) for LLF"),
  make_option("--abc_exclude", type="character", default="AIDS,Depression,HepatitisB,HepatitisC,Migraine,Opioid_use_disorder,Headache,Asthma,Chronic_viral_hepatitis",
              help="Comma-separated disease base names to exclude (without suffix) for ABC"),
  make_option("--llf_min_total", type="double", default=12.9,
              help="Min total cases (incident+prevalent) to display for LLF [default %default]"),
  make_option("--abc_min_total", type="double", default=3,
              help="Min total cases (incident+prevalent) to display for ABC [default %default]"),
  make_option("--llf_title", type="character",
              default="Total Number of Cases per Disease in MGBB-LLF",
              help="Title for LLF plot [default %default]"),
  make_option("--abc_title", type="character",
              default="Total Number of Cases per Disease in MGBB-ABC",
              help="Title for ABC plot [default %default]")
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
      title = title_txt,
      subtitle = subtitle_txt,
      x = "Disease",
      y = "Number of Cases",
      fill = "Case Type"
    ) +
    scale_fill_manual(values = c("Prevalent" = "#1f78b4", "Incident" = "#F8766D")) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      axis.text   = element_text(size = 13),
      axis.title  = element_text(size = 16, face = "bold"),
      plot.title  = element_text(size = 18, face = "bold"),
      legend.title = element_text(size = 14),
      legend.text  = element_text(size = 12)
    )
  ggsave(outfile, p, width = 10, height = 6.5, dpi = 300)
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
# ABC (MGBB-ABC)
# -------------------------------
abc <- fread(opt$abc_phe)
n_abc <- nrow(abc)

exclude_abc <- strsplit(opt$abc_exclude, ",")[[1]] |> trimws()
counts_abc  <- compute_counts(abc, exclude_abc, n_label = n_abc)

# Keep > threshold
counts_abc_fil <- counts_abc %>% filter(total_cases > opt$abc_min_total)

# Save tables
fwrite(counts_abc,      file.path(opt$out_dir, "abc_case_counts_all.csv"))
fwrite(counts_abc_fil,  file.path(opt$out_dir, "abc_case_counts_filtered.csv"))

# Plot Figure 1D
abc_long <- to_long(counts_abc_fil)
plot_counts(
  abc_long,
  title_txt   = sprintf("%s (N=%d)", opt$abc_title, n_abc),
  subtitle_txt= NULL,
  outfile     = file.path(opt$out_dir, "figure1D_abc_disease_counts.png")
)

# Optional summaries
capture.output(summary(abc$COLLECTION_DATE),
               file = file.path(opt$out_dir, "abc_collection_date_summary.txt"))
capture.output(summary(abc$FU_endDate),
               file = file.path(opt$out_dir, "abc_fu_enddate_summary.txt"))
capture.output(summary(abc$FU_years),
               file = file.path(opt$out_dir, "abc_fu_years_summary.txt"))

message("Done. Wrote figures + tables to: ", normalizePath(opt$out_dir))
