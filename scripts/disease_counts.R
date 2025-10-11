#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(argparse); library(data.table); library(dplyr); library(ggplot2)
})

parser <- ArgumentParser()
parser$add_argument("--phe", required=TRUE)
parser$add_argument("--out_prefix", required=TRUE)
args <- parser$parse_args()

phe <- data.table::fread(args$phe)

prev <- grep("_at_collect$", names(phe), value=TRUE)
inc  <- grep("_flag$", names(phe), value=TRUE)
dx <- unique(gsub("(_at_collect|_flag)$", "", c(prev, inc)))

df <- lapply(dx, function(d){
  pcol <- paste0(d, "_at_collect")
  icol <- paste0(d, "_flag")
  pc <- if (pcol %in% names(phe)) sum(phe[[pcol]]==1, na.rm=TRUE) else 0
  ic <- if (icol %in% names(phe)) sum(phe[[icol]]==1, na.rm=TRUE) else 0
  data.frame(disease=d, prevalent=pc, incident=ic)
}) %>% bind_rows()

df_long <- tidyr::pivot_longer(df, cols=c(prevalent, incident), names_to="type", values_to="cases")

p <- ggplot(df_long, aes(x=reorder(disease, cases), y=cases, fill=type))+
  geom_col(position="stack")+
  coord_flip()+
  scale_fill_manual(values=c("prevalent"="#4575b4","incident"="#d73027"))+
  labs(x=NULL, y="Cases", title="Disease counts (prevalent vs incident)")+
  theme_bw(base_size=11)
ggsave(paste0(args$out_prefix,".png"), p, width=8, height=10, dpi=300)
