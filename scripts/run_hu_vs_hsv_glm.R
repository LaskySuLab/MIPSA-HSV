#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(parallel)
  source("R/utils_common.R")
  source("R/utils_qtl.R")
})

opt <- OptionParser(option_list = list(
  make_option("--cohort",       type="character", help="MGBB-LLF or MGBB-ABC"),
  make_option("--phenotype",    type="character", help="Phenotype CSV (cohort cleaned)"),
  make_option("--virus-bin",    type="character", help="TSV from build_binary_matrices.R (HSV)"),
  make_option("--human-bin",    type="character", help="TSV from build_binary_matrices.R (Hu FL)"),
  make_option("--aa",           type="integer",   help="Start index (1-based) of human vars to run"),
  make_option("--bb",           type="integer",   help="End index (1-based) of human vars to run"),
  make_option("--outdir",       type="character", default="results/glm")
)) |> parse_args()

cfg <- cohort_preset(opt$cohort)
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

phe <- fread(opt$phenotype)
phe <- detect_and_unify_id(phe)

vb <- fread(opt$virus_bin)
hb <- fread(opt$human_bin)

vb <- detect_and_unify_id(vb)
hb <- detect_and_unify_id(hb)

# Merge (right joins keep binaries' samples)
DT <- phe %>% right_join(vb, by="Sample_ID") %>% right_join(hb, by="Sample_ID") %>% as.data.table()
DT <- recode_covariates(DT)

hu_vars <- names(hb)[!names(hb) %in% "Sample_ID"]
vi_vars <- names(vb)[!names(vb) %in% "Sample_ID"]

aa <- max(1L, opt$aa)
bb <- min(length(hu_vars), opt$bb)
message("Running human indices: ", aa, ":", bb, " of ", length(hu_vars))

num_cores <- max(1L, detectCores() - 1L)
cl <- makeCluster(num_cores)
clusterExport(cl, list("DT","hu_vars","vi_vars"), envir=environment())
clusterEvalQ(cl, { library(data.table); library(stats) })

run_glm <- function(ii) {
  abs <- hu_vars[ii]
  res_list <- lapply(vi_vars, function(vv) {
    fml <- as.formula(paste0(abs, " ~ ", vv,
                             " + Age_at_collect + Gender + BMI_most_recent + Race_Group + Smoking + Alcohol"))
    fit <- tryCatch(glm(fml, family = "binomial", data = DT), error=function(e) NULL)
    if (is.null(fit)) return(c(NA,NA,NA,NA))
    sm <- summary(fit)
    # coefficient 2 is the virus term
    c(sm$coefficients[2, 1:4])
  })
  out <- as.data.table(do.call(rbind, res_list))
  setnames(out, c("Beta","SE","Z","P"))
  out[, var_id := vi_vars]
  out[, antibody := abs]
  out
}

chunks <- parLapply(cl, aa:bb, run_glm)
stopCluster(cl)

for (i in seq_along(chunks)) {
  abs <- hu_vars[aa + i - 1L]
  fp <- file.path(opt$outdir, sprintf("%s_glm_%s.csv", opt$cohort, abs))
  fwrite(chunks[[i]][order(P)], fp)
}



Rscript scripts/run_hu_vs_hsv_glm.R \
  --cohort MGBB-LLF \
  --phenotype data/MIPSA_Asthma_1290.csv \
  --virus-bin results/binaries/MGBB-LLF_hsv_binary.tsv \
  --human-bin results/binaries/MGBB-LLF_human_fl_binary.tsv \
  --aa 1 --bb 2258
