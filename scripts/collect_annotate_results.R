#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  source("R/utils_common.R")
  source("R/utils_qtl.R")
})

opt <- OptionParser(option_list = list(
  make_option("--cohort",          type="character", help="MGBB-LLF or MGBB-ABC"),
  make_option("--glm-dir",         type="character", help="Directory with per-antibody GLM CSVs"),
  make_option("--virsight-promax", type="character", help="VirSIGHT *Promax_Hits_FOB.csv"),
  make_option("--husight-fl",      type="character", help="HuSIGHT FullLength *Hits_FOB.csv"),
  make_option("--virus-bin",       type="character", help="binary TSV (HSV) used"),
  make_option("--human-bin",       type="character", help="binary TSV (Hu FL) used"),
  make_option("--replicate-of",    type="character", default="", help="(Optional) path to OTHER cohort annotated .rds for replication"),
  make_option("--outdir",          type="character", default="results/annot")
)) |> parse_args()

cfg <- cohort_preset(opt$cohort)
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

# cat all GLM csvs
fls <- list.files(opt$glm_dir, pattern="\\.csv$", full.names=TRUE)
if (!length(fls)) stop("No GLM csvs found in --glm-dir")
tot <- rbindlist(lapply(fls, fread), use.names=TRUE, fill=TRUE)
tot[, direction_of_effect := ifelse(Beta > 0, "positive", "negative")]

# annotations
vi <- fread(opt$virsight_promax)[, .(UniProt_acc, taxon_genus, taxon_species, product)]
hu <- fread(opt$husight_fl)
if (!"var_id" %in% names(hu)) hu[, var_id := paste0("h", .I)]
hu_anno <- hu[, .(antibody = var_id, gene_symbol, UniProt_acc)]

tot <- tot %>%
  left_join(vi, by=c("var_id"="UniProt_acc")) %>%
  left_join(hu_anno, by="antibody")

# FDR
tot[, P.adj := p.adjust(P, method="fdr")]

# Save master
rds_master <- file.path(opt$outdir, sprintf("%s_glm_annotated.rds", opt$cohort))
fwrite(tot, file.path(opt$outdir, sprintf("%s_glm_annotated.csv", opt$cohort)))
saveRDS(tot, rds_master)

# Significant
sig <- tot[P.adj < 0.05]
fwrite(sig, file.path(opt$outdir, sprintf("%s_glm_sig_fdr05.csv", opt$cohort)))

# Replication (if provided)
if (nzchar(opt$replicate_of) && file.exists(opt$replicate_of)) {
  other <- readRDS(opt$replicate_of)
  # Replicate exact antibody–virus–direction tuples
  sig[, key := paste(var_id, antibody, direction_of_effect, sep="__")]
  other[, key := paste(var_id, antibody, direction_of_effect, sep="__")]
  rep_keys <- intersect(sig$key, other$key)
  rep_sig <- sig[key %in% rep_keys]
  fwrite(rep_sig, file.path(opt$outdir, sprintf("%s_rep_sig.csv", opt$cohort)))
}

# Prevalence (per human gene & per virus species)
vb <- fread(opt$virus_bin)
hb <- fread(opt$human_bin)

# collapse by gene_symbol
hb_long <- as.data.table(hb)
hb_long <- hb_long[, !"Sample_ID"]
gene_map <- hu_anno[, .(antibody, gene_symbol, HU_Uni=UniProt_acc)]
setDT(gene_map)
setnames(hb_long, names(hb_long), make.unique(names(hb_long)))
hb_names <- names(hb_long)

# antibody-level prevalence
prev_hu <- data.table(antibody = hb_names,
                      count = colSums(hb_long, na.rm = TRUE))
prev_hu <- prev_hu %>%
  left_join(gene_map, by="antibody") %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol), HU_Uni, gene_symbol),
         prevalence = count * 100 / cfg$n) %>%
  group_by(gene_symbol) %>% summarise(Prevalence = max(prevalence, na.rm=TRUE), .groups="drop") %>%
  arrange(desc(Prevalence))
fwrite(prev_hu, file.path(opt$outdir, sprintf("%s_prevalence_human.csv", opt$cohort)))

# virus species prevalence (max across peptides per species)
vb_long <- as.data.table(vb)[, !"Sample_ID"]
prev_vi <- data.table(var_id = names(vb_long),
                      count = colSums(vb_long, na.rm=TRUE)) %>%
  left_join(vi, by=c("var_id"="UniProt_acc")) %>%
  mutate(prevalence = count * 100 / cfg$n) %>%
  group_by(taxon_species) %>% summarise(Prevalence = max(prevalence, na.rm=TRUE), .groups="drop") %>%
  mutate(taxon_species = species_short(taxon_species))
fwrite(prev_vi, file.path(opt$outdir, sprintf("%s_prevalence_virus_species.csv", opt$cohort)))


# LLF
Rscript scripts/collect_annotate_results.R \
  --cohort MGBB-LLF \
  --glm-dir results/glm \
  --virsight-promax data/IB1007_VirSIGHT_Promax_Hits_Fold-Over-Background.csv \
  --husight-fl     data/IB1007_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv \
  --virus-bin results/binaries/MGBB-LLF_hsv_binary.tsv \
  --human-bin results/binaries/MGBB-LLF_human_fl_binary.tsv

# ABC, referencing LLF file for replication
Rscript scripts/collect_annotate_results.R \
  --cohort MGBB-ABC \
  --glm-dir results/glm \
  --virsight-promax data/IB1021_VirSIGHT_Promax_Hits_Fold-Over-Background.csv \
  --husight-fl     data/IB1021_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv \
  --virus-bin results/binaries/MGBB-ABC_hsv_binary.tsv \
  --human-bin results/binaries/MGBB-ABC_human_fl_binary.tsv \
  --replicate-of results/annot/MGBB-LLF_glm_annotated.rds
