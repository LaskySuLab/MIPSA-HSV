#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse); library(data.table); library(dplyr); library(stringr)
  source("R/utils_common.R"); source("R/utils_qtl.R")
})

option_list <- list(
  make_option("--res_dir",type="character", help="Directory of hu_vs_vi_chunk*.csv"),
  make_option("--vi_anno",type="character", help="VirSIGHT FoB CSV"),
  make_option("--hu_anno",type="character", help="HuSIGHT FullLength FoB CSV"),
  make_option("--rep_sig",type="character", default=NULL, help="rep_hsv.sig"),
  make_option("--out_prefix",type="character", default="results/qtl_hu_vs_vi", help="Output prefix")
)
opt <- parse_args(OptionParser(option_list=option_list))

files <- list.files(opt$res_dir, pattern="^hu_vs_vi_chunk.*\\.csv$", full.names=TRUE)
if (!length(files)) stop("No result chunks found")
res <- rbindlist(lapply(files, fread), fill = TRUE)

vi_anno <- fread(opt$vi_anno)[, .(UniProt_acc, taxon_genus, taxon_species, product)]
hu_anno <- fread(opt$hu_anno)
if (!"var_id" %in% names(hu_anno)) hu_anno[, var_id := paste0("h", seq_len(.N))]

out <- res %>%
  left_join(hu_anno[, .(var_id, gene_symbol, UniProt_acc)], by = c("antibody" = "var_id")) %>%
  left_join(vi_anno, by = c("var_id" = "UniProt_acc")) %>%
  mutate(FDR = bh_fdr(P))

if (!is.null(opt$rep_sig)) {
  rep <- fread(opt$rep_sig)
  rep$qtl <- paste0(rep$taxon_species, "_", rep$gene, "_", rep$antibody)
  out$qtl_key <- paste0(out$taxon_species, "_", out$gene_symbol, "_", out$antibody)
  out$replicated <- out$qtl_key %in% rep$qtl
  out$qtl_key <- NULL
}

fwrite(out, paste0(opt$out_prefix, "_all.csv"))
fwrite(out %>% filter(FDR < 0.05), paste0(opt$out_prefix, "_FDR05.csv"))
message("Wrote: ", opt$out_prefix, "_all.csv and _FDR05.csv")
