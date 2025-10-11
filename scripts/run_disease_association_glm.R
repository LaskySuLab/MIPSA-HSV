#!/usr/bin/env Rscript
# Disease ~ ImmuneMarker (Virus peptide or Autoantigen), aggregated prevalent+incident as in manuscript
suppressPackageStartupMessages({
  library(argparse); library(data.table); library(dplyr); library(tidyr)
  source("R/utils_common.R")
})

parser <- ArgumentParser()
parser$add_argument("--phe", required=TRUE, help="Phenotype with disease columns and covariates")
parser$add_argument("--immune_mat", required=TRUE, help="Binary matrix of markers (rows Sample_ID)")
parser$add_argument("--id_col", default="Sample_ID")
parser$add_argument("--disease_suffixes", nargs="+", default=c("_at_collect","_flag"),
                    help="Columns that will be rowSums() to create *_combine")
parser$add_argument("--min_cases", type="integer", default=10)
parser$add_argument("--out", required=TRUE)
args <- parser$parse_args()

phe <- read_any(args$phe)
setnames(phe, old=args$id_col, new="Sample_ID", skip_absent=TRUE)

# Build combined disease columns
status_cols <- grep("_Status", names(phe), value=TRUE)
diseases <- gsub("_Status","", status_cols)
for (dx in diseases) {
  cols <- paste0(dx, args$disease_suffixes)
  cols <- intersect(cols, names(phe))
  if (length(cols) == 0) next
  phe[, paste0(dx,"_combine") := as.integer(rowSums(.SD, na.rm=TRUE) > 0), .SDcols=cols]
}
dx_cols <- grep("_combine$", names(phe), value=TRUE)

# Filter by case count
dx_mat <- phe[, ..dx_cols]
keep_dx <- names(which(colSums(dx_mat, na.rm=TRUE) >= args$min_cases))
if (length(keep_dx) == 0) stop("No diseases pass min_cases.")

covs <- make_covariates(phe)

M <- read_any(args$immune_mat)
setnames(M, old=names(M)[1], new="Sample_ID")

DD <- merge(covs, M, by="Sample_ID")
out <- vector("list", length(keep_dx))

for (i in seq_along(keep_dx)) {
  dx <- keep_dx[i]
  dat <- merge(DD, phe[, .(Sample_ID, y=get(dx))], by="Sample_ID")
  res_i <- lapply(setdiff(names(M), "Sample_ID"), function(mk) {
    fit <- safe_glm(glm_formula(y="y", x=mk), dat)
    if (is.null(fit)) return(NULL)
    s <- summary(fit)$coefficients
    data.frame(disease=dx, marker=mk, beta=s[mk,"Estimate"], se=s[mk,"Std. Error"], z=s[mk,"z value"], p=s[mk,"Pr(>|z|)"])
  })
  out[[i]] <- rbindlist(res_i, fill=TRUE)
}
res <- rbindlist(out, fill=TRUE)
res[, fdr := p.adjust(p, "BH")]
fwrite(res, args$out)
