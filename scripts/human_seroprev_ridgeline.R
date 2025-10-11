#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(dplyr); library(ggridges); library(ggplot2)
})

option_list <- list(
  make_option("--hu_bin", type="character", help="human_fl_binary.txt"),
  make_option("--hu_anno", type="character", help="HuSIGHT FL FoB CSV"),
  make_option("--out_dir", type="character", default="results", help="Output dir")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

hu_bin <- fread(opt$hu_bin)
hu_anno <- fread(opt$hu_anno)
if (!"var_id" %in% names(hu_anno)) hu_anno[, var_id := paste0("h", seq_len(.N))]

abs <- setdiff(colnames(hu_bin), "Subject_Id")
prev <- colMeans(hu_bin[, ..abs], na.rm = TRUE)
tab  <- data.frame(var_id = abs, seroprev = prev) %>%
  left_join(hu_anno[, .(var_id, gene_symbol)], by="var_id") %>%
  mutate(label = ifelse(is.na(gene_symbol), var_id, gene_symbol))

p <- ggplot(tab, aes(x = seroprev, y = 1)) +
  geom_density_ridges(fill="grey80", color="grey40", scale=3, rel_min_height=0.01) +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(x="Human antibody seroprevalence", y="", title="Human Ab seroprevalence") +
  theme_minimal(base_size=12) + theme(axis.text.y = element_blank())
ggsave(file.path(opt$out_dir, "figure2C_hu_ridgeline.png"), p, width=7, height=3, dpi=300)
