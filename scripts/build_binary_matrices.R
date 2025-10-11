#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(argparse); library(data.table); library(dplyr); library(tidyr)
  source("R/utils_common.R")
})

parser <- ArgumentParser()
parser$add_argument("--virus_in", required=TRUE, help="Path to VirSIGHT matrix (rows=Subject_Id/Sample_ID, cols=UniProt)")
parser$add_argument("--human_in", required=TRUE, help="Path to HuSIGHT FL matrix (rows=Subject_Id/Sample_ID, cols=antibody)")
parser$add_argument("--virus_out", required=TRUE, help="Output path for viral binary matrix")
parser$add_argument("--human_out", required=TRUE, help="Output path for human binary matrix")
parser$add_argument("--id_col", default="Subject_Id")
parser$add_argument("--thresh", type="double", default=1.0)
parser$add_argument("--min_prev", type="double", default=0.01)
args <- parser$parse_args()

virus <- read_any(args$virus_in)
human <- read_any(args$human_in)

virus_bin <- binarize_hits(virus, id_col=args$id_col, thresh=args$thresh, min_prev=args$min_prev)
human_bin <- binarize_hits(human, id_col=args$id_col, thresh=args$thresh, min_prev=args$min_prev)

fwrite(virus_bin, args$virus_out)
fwrite(human_bin, args$human_out)
