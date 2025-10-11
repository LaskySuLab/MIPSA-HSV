#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(argparse); library(data.table); library(dplyr); library(ggplot2); library(ggrepel)
})

parser <- ArgumentParser()
parser$add_argument("--virus_bin", required=TRUE)
parser$add_argument("--virus_annot", required=TRUE)
parser$add_argument("--out_prefix", required=TRUE)
args <- parser$parse_args()

V <- fread(args$virus_bin); setnames(V, old=names(V)[1], new="Sample_ID")
A <- fread(args$virus_annot)

keep <- intersect(names(V), unique(A$UniProt_acc))
prev <- data.frame(UniProt_acc=keep, prev=colMeans(V[, ..keep], na.rm=TRUE))
prev <- prev %>% left_join(A, by="UniProt_acc")

p <- ggplot(prev, aes(x=reorder(UniProt_acc, prev), y=prev, color=taxon_species))+
  geom_point(alpha=0.7)+
  coord_flip()+
  scale_y_continuous(labels=scales::percent)+
  labs(x="Viral peptide (UniProt)", y="Prevalence", color="Species", title="Herpesviridae peptide prevalence")+
  theme_bw(base_size=11)
ggsave(paste0(args$out_prefix,".lollipop.png"), p, width=9, height=10, dpi=300)
