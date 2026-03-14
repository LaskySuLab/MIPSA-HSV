---
title: "Seq_homology_figures"
output: html_document
date: "2026-01-18"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(readr)
library(tidyverse)
library(ggplot2)
library(readxl)
options(bitmapType='cairo')
```


```{r}
# load dataframe 

# load the updated excel file
all_updated_peptides <- read_excel("./MIPSA_pairwise_analysis_llf_updated.xlsx",
                               sheet = "SupTable3.Replicated_FDR",
                               col_names = TRUE)

colnames(all_updated_peptides) <- all_updated_peptides[1,]
all_updated_peptides <- all_updated_peptides[-1, 1:7 ]

```


```{r}
# define plot obs against null density functions 
plotNullvsObs <- function(alignment_obs, alignment_null, title = ""){
  
  plot_df <- bind_rows(
    alignment_obs %>% transmute(score, group = "Observed"),
    alignment_null %>% transmute(score, group = "Null")
  ) %>%
    mutate(group = factor(group, levels = c("Observed", "Null")))
  
  ggplot(plot_df, aes(y = score, x = group)) +
    
    ## distribution
    geom_violin(aes(fill = group),
                trim = FALSE, alpha = 0.4, color = NA) +
    
    labs(y = "Alignment Score", x = "", title = title) +
    theme_minimal(base_size = 14) +
  theme(legend.title = element_blank())
}


```

# define other helper functions
```{r}
extract_water_alignment_info <- function(filepath) {
  lines <- readLines(filepath)
  
  # Extract query and subject IDs
  query_line <- grep("^# 1:", lines, value = TRUE)
  subject_line <- grep("^# 2:", lines, value = TRUE)
  query_id <- sub("^# 1:\\s*", "", query_line)
  subject_id <- sub("^# 2:\\s*", "", subject_line)
  
  # Extract alignment metrics
  score_line <- grep("^# Score:", lines, value = TRUE)
  identity_line <- grep("^# Identity:", lines, value = TRUE)
  similarity_line <- grep("^# Similarity:", lines, value = TRUE)
  gaps_line <- grep("^# Gaps:", lines, value = TRUE)
  # Extract alignment length
  length_line <- grep("^# Length:", lines, value = TRUE)
  alignment_length <- as.numeric(sub("^# Length:\\s*", "", length_line))
  
  score <- as.numeric(sub("^# Score:\\s*", "", score_line))
  
  extract_stat <- function(line) {
    m <- regmatches(line, regexec("(\\d+)/(\\d+) \\(([^%]+)%\\)", line))[[1]]
    c(count = as.numeric(m[2]), total = as.numeric(m[3]), percent = as.numeric(m[4]))
  }
  
  identity <- extract_stat(identity_line)
  similarity <- extract_stat(similarity_line)
  gaps <- extract_stat(gaps_line)
  
  # Extract subject alignment lines robustly
  align_lines <- lines[!grepl("^#", lines) & grepl(paste0("^", subject_id, "\\b"), lines)]
  #print(align_lines)
  if (length(align_lines) == 0) {
    return(data.frame(
      file = basename(filepath),
      query_id = query_id,
      subject_id = subject_id,
      subject_start = NA,
      subject_end = NA,
      score = score,
      alignment_length = alignment_length,
      identity_count = identity["count"],
      identity_percent = identity["percent"],
      similarity_count = similarity["count"],
      similarity_percent = similarity["percent"],
      gaps_count = gaps["count"],
      gaps_percent = gaps["percent"]
    ))
  }
  
  # Extract subject start and end positions properly
  # Get first and last subject alignment lines
  first_line <- align_lines[1]
  last_line <- align_lines[length(align_lines)]
  
  # Extract first numeric value for subject start
  subject_start <- as.numeric(sub("^\\S+\\s+(\\d+).*", "\\1", first_line))
  
  # Extract last numeric value for subject end
  subject_end <- as.numeric(sub(".*\\s(\\d+)\\s*$", "\\1", last_line))
  
  data.frame(
    file = basename(filepath),
    query_id = query_id,
    subject_id = subject_id,
    subject_start = subject_start,
    subject_end = subject_end,
    score = score,
    alignment_length = alignment_length,
    identity_count = identity["count"],
    identity_percent = identity["percent"],
    similarity_count = similarity["count"],
    similarity_percent = similarity["percent"],
    gaps_count = gaps["count"],
    gaps_percent = gaps["percent"]
  )
}

alignment_metrics_summary <- function(alignment_summary){
  summaries <- alignment_summary |>
    summarise(score_mean = mean(score), score_sd = sd(score),
              length_mean = mean(alignment_length), length_sd = sd(alignment_length),
              identity_mean = mean(identity_percent), identity_sd = sd(identity_percent), 
              similarity_mean = mean(similarity_percent), similarity_sd = sd(similarity_percent),
              gaps_mean = mean(gaps_percent), gaps_sd = sd(gaps_percent))
  return(cbind(subject_id = unique(alignment_summary$subject_id),
               as.data.frame(summaries)))
}


```


# CFHR1 vs. CMV
```{r}
# 89 peptides in total
alignment_dir <- "./CFHR1_CMV/water_alignments_assoc/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_CFHR1 <- do.call(rbind, lapply(files, extract_water_alignment_info))

alignment_dir <- "./CFHR1_CMV/water_alignments_randomized_null_assoc/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_CFHR1_null <- do.call(rbind, lapply(files, extract_water_alignment_info))
plotNullvsObs(alignment_summary_assoc_CFHR1, alignment_summary_assoc_CFHR1_null,
              title = "CMV peptide alignments against CFHR1")
```


# DNAJC12 vs. CMV
```{r}
alignment_dir <- "./DNAJC12_CMV/water_alignments_assoc/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_DNAJC12 <- do.call(rbind, lapply(files, extract_water_alignment_info))
# 213 total

alignment_dir <- "./DNAJC12_CMV/water_alignments_randomized_null_assoc/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_DNAJC12_null <- do.call(rbind, lapply(files, extract_water_alignment_info))
alignment_summary_assoc_DNAJC12_null <- alignment_summary_assoc_DNAJC12_null |>
  mutate(query_num = as.numeric(str_match(file, "query_(\\d+)_")[,2])) |>
  filter(query_num >= 1 & query_num <= 213)
# 213000 in total

plotNullvsObs(alignment_summary_assoc_DNAJC12, alignment_summary_assoc_DNAJC12_null,
              title = "CMV peptide alignments against DNAJC12")
```

# LY9 vs. CMV
```{r}
# 3 peptides in total
alignment_dir <- "./LY9_CMV/water_alignments_assoc/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_LY9 <- do.call(rbind, lapply(files, extract_water_alignment_info))

alignment_dir <- "./LY9_CMV/water_alignments_randomized_null_assoc/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_LY9_null <- do.call(rbind, lapply(files, extract_water_alignment_info))
plotNullvsObs(alignment_summary_assoc_LY9, alignment_summary_assoc_LY9_null,
              title = "CMV peptide alignments against LY9")
```


# HEXIM2 vc. CMV
```{r}
# 77 peptides in total
alignment_dir <- "./HEXIM2_CMV/water_alignments_assoc/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_HEXIM2 <- do.call(rbind, lapply(files, extract_water_alignment_info))

alignment_dir <- "./HEXIM2_CMV/water_alignments_randomized_null_assoc/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_HEXIM2_null <- do.call(rbind, lapply(files, extract_water_alignment_info))
plotNullvsObs(alignment_summary_assoc_HEXIM2, alignment_summary_assoc_HEXIM2_null,
              title = "CMV peptide alignments against HEXIM2")
```


# IQCB1 vs. CMV
```{r}
# 205 peptides in total
alignment_dir <- "./IQCB1_CMV/water_alignments_assoc/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_IQCB1 <- do.call(rbind, lapply(files, extract_water_alignment_info))

alignment_dir <- "./IQCB1_CMV/water_alignments_randomized_null_assoc/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_IQCB1_null <- do.call(rbind, lapply(files, extract_water_alignment_info))
plotNullvsObs(alignment_summary_assoc_IQCB1, alignment_summary_assoc_IQCB1_null,
              title = "CMV peptide alignments against IQCB1")
```


# R3HDM2 vc. CMV
```{r}
# 34 peptides in total
alignment_dir <- "./R3HDM2_CMV/water_alignments_assoc/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_R3HDM2<- do.call(rbind, lapply(files, extract_water_alignment_info))

alignment_dir <- "./R3HDM2_CMV/water_alignments_randomized_null_assoc/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_R3HDM2_null <- do.call(rbind, lapply(files, extract_water_alignment_info))
plotNullvsObs(alignment_summary_assoc_R3HDM2, alignment_summary_assoc_R3HDM2_null,
              title = "CMV peptide alignments against R3HDM2")
```

# SLC30A4 vc. CMV
```{r}
# 97 peptides in total
alignment_dir <- "./SLC30A4_CMV/water_alignments_assoc/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_SLC30A4<- do.call(rbind, lapply(files, extract_water_alignment_info))

alignment_dir <- "./SLC30A4_CMV/water_alignments_randomized_null_assoc/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_SLC30A4_null <- do.call(rbind, lapply(files, extract_water_alignment_info))
plotNullvsObs(alignment_summary_assoc_SLC30A4, alignment_summary_assoc_SLC30A4_null,
              title = "CMV peptide alignments against SLC30A4")
```

# ZNF550 vs. CMV
```{r}
alignment_dir <- "./ZNF550_CMV/water_alignments_assoc_2/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_ZNF550<- do.call(rbind, lapply(files, extract_water_alignment_info))

alignment_dir <- "./ZNF550_CMV/water_alignments_randomized_null_assoc_2/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_ZNF550_null <- do.call(rbind, lapply(files, extract_water_alignment_info))
plotNullvsObs(alignment_summary_assoc_ZNF550, alignment_summary_assoc_ZNF550_null,
              title = "CMV peptide alignments against ZNF550")
```

# KCNMB3 vs. EBV
```{r}
# 8 peptides total
alignment_dir <- "./KCNMB3_EBV/water_alignments_assoc/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_KCNMB3 <- do.call(rbind, lapply(files, extract_water_alignment_info))

alignment_dir <- "./KCNMB3_EBV/water_alignments_randomized_null_assoc/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_KCNMB3_null <- do.call(rbind, lapply(files, extract_water_alignment_info))
plotNullvsObs(alignment_summary_assoc_KCNMB3, alignment_summary_assoc_KCNMB3_null,
              title = "EBV peptide alignments against KCNMB3")
```

# PSMB6 vs. EBV
```{r}
# 1 peptide in total
alignment_dir <- "./PSMB6_EBV/water_alignments_assoc/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_PSMB6 <- do.call(rbind, lapply(files, extract_water_alignment_info))

alignment_dir <- "./PSMB6_EBV/water_alignments_randomized_null_assoc/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_PSMB6_null <- do.call(rbind, lapply(files, extract_water_alignment_info))
plotNullvsObs <- function(alignment_obs, alignment_null, title = ""){

  plot_df <- dplyr::bind_rows(
    alignment_obs  %>% dplyr::transmute(score, group = "Observed"),
    alignment_null %>% dplyr::transmute(score, group = "Null")
  ) %>%
    dplyr::mutate(group = factor(group, levels = c("Observed", "Null")))

  ggplot(plot_df, aes(x = group, y = score)) +

    ## violin for both groups (this drives the legend)
    geom_violin(
      aes(fill = group),
      trim = FALSE,
      alpha = 0.4,
      color = NA,
      show.legend = TRUE
    ) +

    ## jitter ONLY for Observed (hide legend so you don't get a 2nd one)
    geom_jitter(
    data = dplyr::filter(plot_df, group == "Observed"),
    color = "#F8766D",
    width = 0.1,
    size = 2,
    alpha = 0.7
  )+

    scale_fill_manual(
      values = c("Observed" = "#F8766D", "Null" = "#00BFC4"),
      drop = FALSE
    ) +

    labs(y = "Alignment Score", x = "", title = title, fill = NULL) +
    theme_minimal(base_size = 14)+
    theme(legend.title = element_blank())
}
plotNullvsObs(alignment_summary_assoc_PSMB6, alignment_summary_assoc_PSMB6_null,
              title = "EBV peptide alignments against PSMB6")

```

# PHLDA1 vs. HSV1
```{r}
# 232 peptides in total
alignment_dir <- "./PHLDA1_HSV1/water_alignments_assoc/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_phlda1 <- do.call(rbind, lapply(files, extract_water_alignment_info))

alignment_dir <- "./PHLDA1_HSV1/water_alignments_randomized_null_assoc/"
files <- list.files(alignment_dir, pattern = "\\.txt$", full.names = TRUE)
alignment_summary_assoc_phlda1_null <- do.call(rbind, lapply(files, extract_water_alignment_info))
plotNullvsObs(alignment_summary_assoc_phlda1, alignment_summary_assoc_phlda1_null,
              title = "HSV1 peptide alignments against PHLDA1")
```


# create a summary variable (proportion of peptides reaching the top 5% emprical null distribution)
```{r}
prop_in_null <- function(null, obs, name){
  pct <- ecdf(null$score)(obs$score)
  ct <- sum(pct > 0.95)
  prop <- sum(pct > 0.95)/length(pct)
  num_sig_association = nrow(obs)
  if(ct == 0){
    range_quantile_tophits <- NA
  }else{
    range_quantile_tophits <- range(pct[pct > 0.95])
    range_quantile_tophits <- paste0(round(range_quantile_tophits[1], 3), "-", round(range_quantile_tophits[2], 3))
  }
  prop <- cbind(alignment = name, 
                count_95_null = ct,
                num_sig_association = num_sig_association,
                prop_95_null = prop,
                range_quantile_tophits = range_quantile_tophits)
  return(as.data.frame(prop))
}

prop_total <- rbind(
prop_in_null( alignment_summary_assoc_CFHR1_null, alignment_summary_assoc_CFHR1, "CMV_CFHR1"),
prop_in_null( alignment_summary_assoc_DNAJC12_null, alignment_summary_assoc_DNAJC12, "CMV_DNAJC12"),
prop_in_null( alignment_summary_assoc_GABRE1_null, alignment_summary_assoc_GABRE1, "CMV_IQCB1"),
prop_in_null( alignment_summary_assoc_HEXIM2_null, alignment_summary_assoc_HEXIM2, "CMV_HEXIM2"),
prop_in_null(alignment_summary_assoc_LY9_null, alignment_summary_assoc_LY9, "CMV_LY9"),
prop_in_null( alignment_summary_assoc_R3HDM2_null, alignment_summary_assoc_R3HDM2, "CMV_R3HDM2"),
prop_in_null( alignment_summary_assoc_SLC30A4_null, alignment_summary_assoc_SLC30A4, "CMV_SLC30A4"),
prop_in_null( alignment_summary_assoc_ZNF550_null, alignment_summary_assoc_ZNF550, "CMV_ZNF550"),
prop_in_null( alignment_summary_assoc_KCNMB3_null, alignment_summary_assoc_KCNMB3, "EBV_KCNMB3"),
prop_in_null( alignment_summary_assoc_P3H4_null, alignment_summary_assoc_PSMB6, "EBV_PSMB6"),
prop_in_null( alignment_summary_assoc_phlda1_null, alignment_summary_assoc_phlda1, "HSV1_PHLDA1")
)

prop_total
```


```{r}

alignment_summaries_all_peptides <- do.call(rbind, lapply(list(                       alignment_summary_assoc_CFHR1,
                                                               alignment_summary_assoc_DNAJC12,
                                                               alignment_summary_assoc_IQCB1,
                                                               alignment_summary_assoc_HEXIM2,
                                                            alignment_summary_assoc_R3HDM2,
                                                            alignment_summary_assoc_LY9,
                                                               alignment_summary_assoc_SLC30A4,
                                                               alignment_summary_assoc_ZNF550,
                                                               alignment_summary_assoc_KCNMB3,
                                                            alignment_summary_assoc_PSMB6,                                   
                                                               alignment_summary_assoc_phlda1), 
                                                          alignment_metrics_summary))

alignment_summaries_all_peptides <- cbind(prop_total, alignment_summaries_all_peptides)
alignment_summaries_all_peptides <- alignment_summaries_all_peptides |>
  select(-c("subject_id"))

alignment_summaries_all_peptides
write.csv(alignment_summaries_all_peptides, "./alignment_summaries_all_peptides_UPDATED.csv")
alignment_summaries_all_peptides
```

# select top hits 
```{r}
# CMV_ZNF550
CMV_ZNF550_top_hits <- alignment_summary_assoc_ZNF550[order(alignment_summary_assoc_ZNF550$score, decreasing = TRUE),]
CMV_ZNF550_top_hits <- CMV_ZNF550_top_hits[c(1:119), ]


CMV_ZNF550_top_hits |>
  ggplot(aes(y = query_id)) +
  geom_segment(
    aes(x = subject_start, xend = subject_end, yend = query_id),
    linewidth = 0.8,
  ) +
  labs(
    x = "Subject (ZNF550) position",
    y = "CMV peptides"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_blank())
```

```{r}
# extract files that are top hits in the water alignments, ZNF 550 vs. CMV
files <- CMV_ZNF550_top_hits$file
prefix <- sub("_vs.*$", "", files)
fasta_names <- paste0(prefix, ".fasta")
write.table(fasta_names, "~/MIPSA_sequence_homology/ZNF550_CMV/fasta_names.txt")

fastadir <- "~/MIPSA_sequence_homology/ZNF550_CMV/water_alignments_assoc_2"
paths <- file.path(fastadir, unique(fasta_names))

out <- "~/MIPSA_sequence_homology/ZNF550_CMV/combined_queries.fasta"
cat(file = out)                  # empty / create
for (p in paths[file.exists(paths)]) {
  print(p)
  cat(readLines(p), sep = "\n", file = out, append = TRUE)
  cat("\n", file = out, append = TRUE)
}

```



```{r}
# extract files that are top hits in the water alignments, PHLDA1 vs. HSV1
PHLDA1_HSV1_top_hits <- alignment_summary_assoc_phlda1[order(alignment_summary_assoc_phlda1$score, decreasing = TRUE),]
PHLDA1_HSV1_top_hits <- PHLDA1_HSV1_top_hits[c(1:38), ]


PHLDA1_HSV1_top_hits |>
  ggplot(aes(y = query_id)) +
  geom_segment(
    aes(x = subject_start, xend = subject_end, yend = query_id),
    linewidth = 0.8,
  ) +
  labs(
    x = "Subject (PHLDA1) position",
    y = "HSV1 peptides"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_blank())
```
```{r}
# extract files that are top hits in the water alignments, ZNF 550 vs. CMV
files <- PHLDA1_HSV1_top_hits$file
prefix <- sub("_vs.*$", "", files)
fasta_names <- paste0(prefix, ".fasta")
write.table(fasta_names, "~/MIPSA_sequence_homology/PHLDA1_HSV1/fasta_names.txt")

fastadir <- "~/MIPSA_sequence_homology/PHLDA1_HSV1/water_alignments_assoc"
paths <- file.path(fastadir, unique(fasta_names))

out <- "~/MIPSA_sequence_homology/PHLDA1_HSV1/combined_queries.fasta"
cat(file = out)                  # empty / create
for (p in paths[file.exists(paths)]) {
  print(p)
  cat(readLines(p), sep = "\n", file = out, append = TRUE)
  cat("\n", file = out, append = TRUE)
}
```

