#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(argparse); library(data.table); library(dplyr); library(ggplot2); library(ggrepel)
})

parser <- ArgumentParser()
parser$add_argument("--replicated_pairs", required=TRUE, help="Output of collect_annotate_results.R")
parser$add_argument("--out_prefix", required=TRUE)
args <- parser$parse_args()

X <- data.table::fread(args$replicated_pairs)

# DIAMOND: count per autoantigen by species, color by best P
d1 <- X[, .(n_peptides=.N, bestP=min(p_r, na.rm=TRUE)), by=.(taxon_species, hu_gene)]
p1 <- ggplot(d1, aes(x=reorder(hu_gene, -n_peptides), y=n_peptides, size=n_peptides, color=-log10(bestP)))+
  geom_point(shape=18)+
  facet_wrap(~taxon_species, scales="free_x")+
  scale_color_viridis_c(name="-log10(P)")+
  labs(x="Autoantigen", y="# replicated peptides", title="Replicated viral→autoantigen associations")+
  theme_bw(base_size=11)+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
ggsave(paste0(args$out_prefix,".diamond.png"), p1, width=12, height=6, dpi=300)
