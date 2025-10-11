#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(dplyr); library(purrr); library(stringr)
  source("R/utils_common.R"); source("R/utils_qtl.R")
})

option_list <- list(
  make_option("--phe",      type="character", help="Phenotype CSV (MIPSA_Asthma_1290.csv)"),
  make_option("--vi_bin",   type="character", help="virus_proteins_binary.txt"),
  make_option("--hu_bin",   type="character", help="human_fl_binary.txt"),
  make_option("--vi_anno",  type="character", help="VirSIGHT FoB CSV (for species/product)"),
  make_option("--aa",       type="integer", default=1, help="chunk id (1-based)"),
  make_option("--bb",       type="integer", default=1, help="total chunks"),
  make_option("--species_list", type="character", default=NULL,
              help="Comma-separated virus species to include (default: all in vi_anno)"),
  make_option("--out_dir",  type="character", default="results/Res_hu_vs_vi", help="Output dir"),
  make_option("--rep_sig",  type="character", default=NULL, help="Optional rep_hsv.sig for filtering/labels")
)
opt <- parse_args(OptionParser(option_list=option_list))
ensure_dir(opt$out_dir)

phe <- fread(opt$phe); phe[, Subject_Id := as.character(Subject_Id)]
phe <- build_disease_combine(phe)
phe <- recode_covariates(phe)

vi_bin <- fread(opt$vi_bin)
hu_bin <- fread(opt$hu_bin)

stopif_missing_cols(hu_bin, "Subject_Id"); stopif_missing_cols(vi_bin, "Subject_Id")
stopif_missing_cols(phe, "Sample_ID")

vi_anno <- fread(opt$vi_anno)
if (!"var_id" %in% names(hu_bin)) {
  # columns are already h#### as in build step
}

# select virus set (species filter)
if (!is.null(opt$species_list)) {
  keep_species <- str_split(opt$species_list, ",", simplify = TRUE)
  vi_keep <- vi_anno %>% filter(taxon_species %in% keep_species) %>% pull(UniProt_acc) %>% unique()
} else {
  vi_keep <- unique(vi_anno$UniProt_acc)
}
vi_keep <- intersect(vi_keep, colnames(vi_bin))

# merge phe + hu + virus
dt <- phe %>%
  right_join(hu_bin, by = c("Sample_ID" = "Subject_Id")) %>%
  right_join(vi_bin[, c("Subject_Id", vi_keep), with=FALSE], by = c("Sample_ID" = "Subject_Id")) %>%
  as.data.table()

genes <- intersect(colnames(hu_bin)[colnames(hu_bin)!="Subject_Id"], colnames(dt))
virus_p <- vi_keep

# chunk virus peptides
chunks <- split(virus_p, cut(seq_along(virus_p), breaks = opt$bb, labels = FALSE))
virus_chunk <- chunks[[opt$aa]]
if (length(virus_chunk) == 0) stop("Empty virus chunk for --aa=", opt$aa)

covar <- c("Age_at_collect", "Gender", "BMI_most_recent", "Race_Group", "Smoking", "Alcohol")
stopif_missing_cols(dt, c("Sample_ID", covar))

res_list <- list()
for (ab in genes) {
  if (length(unique(dt[[ab]])) < 2) next
  for (vp in virus_chunk) {
    if (length(unique(dt[[vp]])) < 2) next
    fml <- as.formula(paste0(ab, " ~ ", vp, " + ", paste(covar, collapse = " + ")))
    s <- glm_safe(fml, dt)
    if (!is.null(s)) {
      res_list[[length(res_list)+1]] <- cbind(
        data.frame(antibody = ab, var_id = vp, stringsAsFactors = FALSE), s
      )
    }
  }
}

res <- bind_rows(res_list)
if (nrow(res)) {
  # annotate
  res <- res %>%
    left_join(vi_anno[, .(UniProt_acc, taxon_genus, taxon_species, product)], by = c("var_id" = "UniProt_acc")) %>%
    mutate(FDR = bh_fdr(P))
}
outfile <- file.path(opt$out_dir, sprintf("hu_vs_vi_chunk%02dof%02d.csv", opt$aa, opt$bb))
fwrite(res, outfile)
message("Wrote: ", outfile)
