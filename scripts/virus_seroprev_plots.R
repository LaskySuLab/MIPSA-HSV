#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(dplyr); library(ggplot2); library(ggridges); library(tidyr)
})

option_list <- list(
  make_option("--vi_bin", type="character", help="virus_proteins_binary.txt"),
  make_option("--vi_anno", type="character", help="VirSIGHT FoB CSV"),
  make_option("--out_dir", type="character", default="results", help="Output dir")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

vi_bin <- fread(opt$vi_bin)
anno   <- fread(opt$vi_anno)

proteins <- setdiff(colnames(vi_bin), "Subject_Id")
prev <- colMeans(vi_bin[, ..proteins], na.rm = TRUE)
tab  <- data.frame(UniProt_acc = proteins, seroprev = prev) %>%
  left_join(anno[, .(UniProt_acc, taxon_species, product)], by="UniProt_acc")

# Lollipop by species
sp_prev <- tab %>% group_by(taxon_species) %>%
  summarise(mean_prev = mean(seroprev), .groups="drop") %>%
  arrange(desc(mean_prev)) %>%
  mutate(taxon_species = factor(taxon_species, levels = rev(taxon_species)))

p1 <- ggplot(sp_prev, aes(taxon_species, mean_prev)) +
  geom_segment(aes(xend=taxon_species, y=0, yend=mean_prev), color="grey70") +
  geom_point(size=3, color="firebrick") +
  coord_flip() + scale_y_continuous(labels = scales::percent_format()) +
  labs(x="", y="Mean seroprevalence", title="Virus species seroprevalence (lollipop)") +
  theme_minimal(base_size=12)
ggsave(file.path(opt$out_dir, "figure2A_species_lollipop.png"), p1, width=6, height=7, dpi=300)

# Ridgeline of per-protein seroprevalence by species
p2 <- ggplot(tab, aes(x = seroprev, y = taxon_species)) +
  geom_density_ridges(scale=3, rel_min_height = 0.01) +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(x="Protein seroprevalence", y="", title="Per-protein seroprevalence by species") +
  theme_minimal(base_size=12)
ggsave(file.path(opt$out_dir, "figure2B_species_ridgeline.png"), p2, width=7, height=7, dpi=300)
