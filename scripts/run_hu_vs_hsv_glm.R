#!/usr/bin/env Rscript
# Core pairwise Hu ~ Virus tests (manuscript primary analysis), plus optional Virus ~ Hu.
suppressPackageStartupMessages({
  library(argparse); library(data.table); library(dplyr); library(tidyr); library(purrr)
  source("R/utils_common.R")
})

parser <- ArgumentParser()
parser$add_argument("--cohort", choices=c("discovery","replication"), required=TRUE)
parser$add_argument("--phe", required=TRUE, help="Phenotype/covariate file with Sample_ID and covariates")
parser$add_argument("--virus_bin", required=TRUE, help="Binary viral matrix (Subject_Id/Sample_ID rows)")
parser$add_argument("--human_bin", required=TRUE, help="Binary human matrix (Subject_Id/Sample_ID rows)")
parser$add_argument("--virus_annot", required=TRUE, help="Viral annotation table with columns: UniProt_acc, taxon_species, taxon_genus, product")
parser$add_argument("--human_annot", required=TRUE, help="Human annotation with columns: antibody, UniProt_acc, gene_symbol")
parser$add_argument("--species_keep", nargs="+", default=c("Human alphaherpesvirus 1","Human alphaherpesvirus 2","Human gammaherpesvirus 4","Human betaherpesvirus 5"))
parser$add_argument("--min_prev", type="double", default=0.01, help="Minimum prevalence for both marker types in this cohort")
parser$add_argument("--reverse_test", action="store_true", help="Also run Virus ~ Hu (optional)")
parser$add_argument("--out_prefix", required=TRUE)
args <- parser$parse_args()

# Read inputs
phe  <- read_any(args$phe)
covs <- make_covariates(phe)

V <- read_any(args$virus_bin)
H <- read_any(args$human_bin)

# Harmonize IDs
setnames(V, old=names(V)[1], new="Sample_ID")
setnames(H, old=names(H)[1], new="Sample_ID")

DD <- Reduce(function(x,y) merge(x,y,by="Sample_ID"), list(covs, V, H))
if (nrow(DD) == 0) stop("No ID overlap among phe/virus/human matrices.")

# Marker lists with prevalence filter
virus_cols <- setdiff(names(V), "Sample_ID")
human_cols <- setdiff(names(H), "Sample_ID")
prevV <- colMeans(DD[, ..virus_cols], na.rm=TRUE); keepV <- names(prevV[prevV >= args$min_prev])
prevH <- colMeans(DD[, ..human_cols], na.rm=TRUE); keepH <- names(prevH[prevH >= args$min_prev])

if (length(keepV)==0 || length(keepH)==0) stop("No markers pass prevalence filter.")

# Annotations
va <- read_any(args$virus_annot) %>% as.data.table()
ha <- read_any(args$human_annot) %>% as.data.table()

# Restrict viral markers to requested species
va_keep <- va[taxon_species %in% args$species_keep]
keepV   <- intersect(keepV, unique(va_keep$UniProt_acc))

message("Testing N(autoantigens)=", length(keepH), " x N(viral peptides)=", length(keepV))

# ---------- Hu ~ Virus (primary) ----------
run_one <- function(h, v) {
  fml <- glm_formula(y=h, x=v)
  fit <- safe_glm(fml, DD)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)$coefficients
  c(beta=s[v,"Estimate"], se=s[v,"Std. Error"], z=s[v,"z value"], p=s[v,"Pr(>|z|)"])
}

grid <- CJ(var_h=keepH, var_v=keepV)
res_list <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  res <- run_one(grid$var_h[i], grid$var_v[i])
  if (is.null(res)) next
  res_list[[i]] <- data.frame(var_id_h=grid$var_h[i], var_id_v=grid$var_v[i], t(res))
}
res_hu_vs_virus <- rbindlist(res_list, fill=TRUE)

# Join annotations
res_hu_vs_virus <- res_hu_vs_virus %>%
  left_join(ha %>% select(antibody, UniProt_acc, gene_symbol) %>% rename(var_id_h=antibody, hu_UniProt=UniProt_acc, hu_gene=gene_symbol), by="var_id_h") %>%
  left_join(va_keep %>% select(UniProt_acc, taxon_genus, taxon_species, product) %>% rename(var_id_v=UniProt_acc), by="var_id_v") %>%
  as.data.table()
res_hu_vs_virus[, fdr := bh(p)]

fwrite(res_hu_vs_virus, sprintf("%s.%s.hu_vs_virus.csv.gz", args$out_prefix, args$cohort))

# ---------- Optional Virus ~ Hu ----------
if (isTRUE(args$reverse_test)) {
  run_rev <- function(v, h) {
    fml <- glm_formula(y=v, x=h)
    fit <- safe_glm(fml, DD)
    if (is.null(fit)) return(NULL)
    s <- summary(fit)$coefficients
    c(beta=s[h,"Estimate"], se=s[h,"Std. Error"], z=s[h,"z value"], p=s[h,"Pr(>|z|)"])
  }
  res_list2 <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    res <- run_rev(grid$var_v[i], grid$var_h[i])
    if (is.null(res)) next
    res_list2[[i]] <- data.frame(var_id_v=grid$var_v[i], var_id_h=grid$var_h[i], t(res))
  }
  res_virus_vs_hu <- rbindlist(res_list2, fill=TRUE) %>%
    left_join(va_keep %>% select(UniProt_acc, taxon_genus, taxon_species, product) %>% rename(var_id_v=UniProt_acc), by="var_id_v") %>%
    left_join(ha %>% select(antibody, UniProt_acc, gene_symbol) %>% rename(var_id_h=antibody, hu_UniProt=UniProt_acc, hu_gene=gene_symbol), by="var_id_h") %>%
    as.data.table()
  res_virus_vs_hu[, fdr := bh(p)]
  fwrite(res_virus_vs_hu, sprintf("%s.%s.virus_vs_hu.csv.gz", args$out_prefix, args$cohort))
}
