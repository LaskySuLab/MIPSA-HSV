#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(dplyr); library(tidyr); library(ggplot2); library(stringr)
})

option_list <- list(
  make_option("--phe", type="character", help="MIPSA_Asthma_1290.csv"),
  make_option("--out_dir", type="character", default="results", help="Output dir")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

phe <- fread(opt$phe)
status_vars  <- grep("_Status$", colnames(phe), value = TRUE)
collect_vars <- grep("_at_collect$", colnames(phe), value = TRUE)
disease <- gsub("_Status$", "", status_vars)

for (prefix in disease) {
  cols <- c(paste0(prefix,"_at_collect"), paste0(prefix,"_flag"))
  cols <- cols[cols %in% names(phe)]
  if (!length(cols)) next
  phe[, (paste0(prefix,"_combine")) := as.integer(rowSums(.SD, na.rm = TRUE) > 0), .SDcols = cols]
}

dx <- grep("_combine$", colnames(phe), value = TRUE)
counts <- data.frame(
  disease = gsub("_combine$", "", dx),
  prevalent_cases = sapply(dx, function(x) sum(phe[[x]]==1, na.rm = TRUE))
)

p <- counts %>%
  arrange(desc(prevalent_cases)) %>%
  head(40) %>%
  mutate(disease = factor(disease, levels = rev(disease))) %>%
  ggplot(aes(disease, prevalent_cases)) +
  geom_col(fill = "steelblue") +
  coord_flip() + labs(x="", y="Cases", title="Top disease counts (combined prevalent+incident)") +
  theme_minimal(base_size = 12)

ggsave(file.path(opt$out_dir, "figure1_disease_counts.png"), p, width=6, height=8, dpi=300)
