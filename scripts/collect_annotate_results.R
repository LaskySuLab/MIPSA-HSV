#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(argparse); library(data.table); library(dplyr); library(tidyr); library(stringr)
  source("R/utils_common.R")
})

parser <- ArgumentParser()
parser$add_argument("--disc_hu_vs_virus", required=TRUE, help="Discovery Hu~Virus results")
parser$add_argument("--repl_hu_vs_virus", required=TRUE, help="Replication Hu~Virus results")
parser$add_argument("--virus_annot", required=TRUE)
parser$add_argument("--human_annot", required=TRUE)
parser$add_argument("--out", required=TRUE)
args <- parser$parse_args()

disc <- read_any(args$disc_hu_vs_virus)
repl <- read_any(args$repl_hu_vs_virus)

# strict ID matching and direction consistency
disc_sig <- disc %>% filter(fdr < 0.05)
rep_join <- repl %>%
  select(var_id_h, var_id_v, beta_r=beta, p_r=p) %>%
  inner_join(disc_sig %>% select(var_id_h, var_id_v, beta_d=beta, p_d=p, fdr), by=c("var_id_h","var_id_v")) %>%
  mutate(direction_match = sign(beta_r) == sign(beta_d),
         replicated = direction_match & p_r < 0.05)

va <- read_any(args$virus_annot); ha <- read_any(args$human_annot)
final <- rep_join %>%
  left_join(ha %>% select(antibody, UniProt_acc, gene_symbol) %>% rename(var_id_h=antibody, hu_UniProt=UniProt_acc, hu_gene=gene_symbol), by="var_id_h") %>%
  left_join(va %>% select(UniProt_acc, taxon_genus, taxon_species, product) %>% rename(var_id_v=UniProt_acc), by="var_id_v")

fwrite(final, args$out)
