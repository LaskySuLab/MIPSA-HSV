#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(argparse); library(data.table); library(dplyr); library(igraph); library(ggplot2); library(ggrepel); library(scales)})
parser <- ArgumentParser(); parser$add_argument("--human_dx", required=TRUE); parser$add_argument("--virus_dx", required=TRUE); parser$add_argument("--rep_qtl", required=TRUE); parser$add_argument("--outdir", required=TRUE)
args <- parser$parse_args(); dir.create(args$outdir, showWarnings=FALSE, recursive=TRUE)
human_dx <- fread(args$human_dx); proid_dx <- fread(args$virus_dx); human_qtl <- fread(args$rep_qtl)
human_dx$gene_symbol  <- ifelse(is.na(human_dx$gene_symbol),  human_dx$UniProt_acc,  human_dx$gene_symbol)
proid_dx$gene_symbol  <- ifelse(is.na(proid_dx$gene_symbol),  proid_dx$UniProt_acc,  proid_dx$gene_symbol)
human_qtl$gene_symbol <- ifelse(is.na(human_qtl$gene_symbol), human_qtl$UniProt_acc, human_qtl$gene_symbol)
human_dx$disease <- gsub("_", " ", human_dx$disease); proid_dx$disease <- gsub("_", " ", proid_dx$disease)
source_network <- function(virus_focus, seed=123, file="network.png", title_prefix=NULL){
  tq <- human_qtl %>% filter(taxon_species==virus_focus); hsub <- human_dx %>% filter(var_id %in% tq$antibody); vsub <- proid_dx %>% filter(var_id %in% tq$var_id)
  edges_vh <- tq %>% filter(!is.na(P), P<0.05) %>% transmute(from=var_id,to=gene_symbol,logp=-log10(P),etype="Virus→Auto-Ab")
  edges_hd <- hsub %>% filter(!is.na(P), P<0.05) %>% transmute(from=gene_symbol,to=disease,logp=-log10(P),etype="Auto-Ab→Disease")
  edges_vd <- vsub %>% filter(!is.na(P), P<0.05) %>% transmute(from=var_id,to=disease,logp=-log10(P),etype="Virus→Disease")
  edges_all <- dplyr::bind_rows(edges_vh,edges_hd,edges_vd); if (!nrow(edges_all)) return(invisible(NULL))
  node_names <- unique(c(edges_all$from, edges_all$to))
  nodes <- data.frame(name=node_names) %>% mutate(ntype=dplyr::case_when(name %in% tq$var_id ~ "Virus", name %in% tq$gene_symbol ~ "Auto-Ab", name %in% unique(c(hsub$disease,vsub$disease)) ~ "Disease", TRUE ~ "Other"))
  g <- graph_from_data_frame(edges_all, vertices=nodes, directed=TRUE); w <- edges_all$logp; w[!is.finite(w) | w<=0] <- 1e-3; set.seed(seed); lay <- layout_with_fr(g, weights=w)
  layout_df <- data.frame(x=lay[,1], y=lay[,2], name=V(g)$name); nodes <- left_join(nodes, layout_df, by="name")
  p <- ggplot() + geom_segment(data=edges_all %>% left_join(nodes, by=c("from"="name")) %>% rename(x_from=x,y_from=y) %>% left_join(nodes, by=c("to"="name")) %>% rename(x_to=x,y_to=y),
      aes(x=x_from,y=y_from,xend=x_to,yend=y_to,color=etype,linewidth=logp), alpha=0.5) +
    scale_color_manual(values=c("Virus→Disease"="palegreen1","Virus→Auto-Ab"="royalblue2","Auto-Ab→Disease"="red"), name="Edge Type") +
    scale_linewidth_continuous(range=c(1,3), name="-log10(p)") +
    geom_point(data=nodes, aes(x=x,y=y,fill=ntype), size=3, shape=21, color="grey20") +
    scale_fill_manual(values=c("Auto-Ab"="skyblue","Virus"="green","Disease"="pink","Other"="grey70"), name="Node Type") +
    geom_text(data=nodes %>% filter(ntype=="Auto-Ab"), aes(x=x,y=y,label=name), size=3.6, fontface="bold") +
    theme_void() + theme(legend.position="bottom") + labs(title=paste0(title_prefix," Network: ", gsub("^Human ","",virus_focus)))
  ggsave(file.path(args$outdir, file), p, width=10, height=10)
}
source_network("Human alphaherpesvirus 1", seed=123, file="network_plot_HSV1.png", title_prefix="HSV-1")
source_network("Human alphaherpesvirus 2", seed=123, file="network_plot_HSV2.png", title_prefix="HSV-2")
source_network("Human betaherpesvirus 5", seed=201, file="network_plot_CMV.png",  title_prefix="CMV")
source_network("Human gammaherpesvirus 4", seed=99,  file="network_plot_EBV.png",  title_prefix="EBV")
