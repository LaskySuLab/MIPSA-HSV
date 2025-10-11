#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(argparse); library(data.table); library(dplyr); library(ggplot2); library(ggridges)
})

parser <- ArgumentParser()
parser$add_argument("--human_bin", required=TRUE)
parser$add_argument("--human_annot", required=TRUE)
parser$add_argument("--out_prefix", required=TRUE)
args <- parser$parse_args()

H <- fread(args$human_bin); setnames(H, old=names(H)[1], new="Sample_ID")
A <- fread(args$human_annot)

keep <- intersect(names(H), unique(A$antibody))
prev <- data.frame(antibody=keep, prev=colMeans(H[, ..keep], na.rm=TRUE))
prev <- prev %>% left_join(A, by="antibody")

p <- ggplot(prev, aes(x=prev, y=..density..))+
  geom_histogram(bins=60, alpha=0.8)+
  scale_x_continuous(labels=scales::percent)+
  labs(x="Autoantigen prevalence", y="Density", title="HuSIGHT full-length autoantigen prevalence")+
  theme_bw(base_size=11)
ggsave(paste0(args$out_prefix,".hist.png"), p, width=7, height=5, dpi=300)
