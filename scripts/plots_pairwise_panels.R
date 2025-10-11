#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(viridis)
  library(ComplexHeatmap)
  library(circlize)
  source("R/utils_common.R")
  source("R/utils_qtl.R")
})

opt <- OptionParser(option_list = list(
  make_option("--cohort",      type="character", help="MGBB-LLF or MGBB-ABC"),
  make_option("--annot-rds",   type="character", help="collect_annotate_results.R output .rds (annotated)"),
  make_option("--prev-human",  type="character", help="human prevalence csv"),
  make_option("--prev-virus",  type="character", help="virus species prevalence csv"),
  make_option("--other-annot-rds", type="character", default="", help="(optional) other cohort annotated rds for matched ordering"),
  make_option("--outdir",      type="character", default="results/plots_qtl")
)) |> parse_args()

cfg <- cohort_preset(opt$cohort)
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)

dat <- readRDS(opt$annot_rds)
dat[, taxon_species := species_short(taxon_species)]
dat[, log10P := -log10(ifelse(is.na(P.adj), P, P.adj))]

# ---- Diamond (bubble) plot per cohort ----
# Size = association count per (species, gene), color = max -log10P
bub <- dat %>%
  group_by(taxon_species, gene_symbol) %>%
  summarise(Association_Count = n(),
            Max_log10_P = max(log10P, na.rm=TRUE),
            .groups="drop") %>%
  left_join(read.csv(opt$prev_human), by="gene_symbol")
bub$taxon_species <- factor(bub$taxon_species, levels=c("HSV-1","HSV-2","VZV","CMV","HHV-6A","HHV-6B","HHV-7","EBV","HHV-8"))

bub$gene_label <- sprintf("%s (%.2f)", bub$gene_symbol, bub$Prevalence)
bub <- bub %>% arrange(desc(Max_log10_P)) %>%
  mutate(gene_label = factor(gene_label, levels = unique(gene_label)))

p_dia <- ggplot(bub, aes(x = gene_label, y = taxon_species)) +
  geom_point(aes(size = pmin(Association_Count, 100), color = pmin(Max_log10_P, ifelse(cfg$title=="MGBB-LLF", 80, 12))),
             shape=18, alpha=0.8) +
  scale_size_continuous(range=c(2,10), breaks=c(1,3,5,10,30), labels=c("1","3","5","10","30+")) +
  scale_color_viridis_c(option="viridis") +
  labs(
    title = sprintf("Significant Auto-antibodies associated with Herpesviruses in %s", cfg$title),
    x = "Gene Symbols (prevalence %)", y = "Herpes viruses",
    size = "Association Count", color = "-log10(P)"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=9),
        plot.title = element_text(hjust=0.5, size=15))
ggsave(file.path(opt$outdir, sprintf("diamond_%s.png", cfg$title)), p_dia, width=12, height=5, dpi=300)

# ---- Manhattan-like panel across species (collapse per antibody x product min P) ----
dat[, product1 := paste(antibody, product, sep="_")]
coll <- dat %>% group_by(product1) %>% slice_min(order_by = P, with_ties = FALSE) %>% ungroup()
coll$taxon_species <- factor(coll$taxon_species, levels=c("HSV-1","HSV-2","VZV","CMV","HHV-6A","HHV-6B","HHV-7","EBV","HHV-8"))
p_thr <- if (cfg$title=="MGBB-LLF")  -log10(max(dat[P.adj<0.05, P, na.rm=TRUE])) else -log10(max(dat[p.adjust(P, "fdr")<0.05, P, na.rm=TRUE]))

p_manh <- ggplot(coll, aes(x=antibody, y = -log10(P))) +
  facet_grid(. ~ taxon_species, scales="free_x") +
  geom_point(aes(color=Beta, shape=direction_of_effect), size=2.2, alpha=1) +
  scale_shape_manual(values=c("negative"=25,"positive"=17)) +
  scale_color_gradient2(low="#2b83ba", mid="white", high="#d7191c", midpoint=0, limits=c(-5,5), oob=scales::squish) +
  geom_hline(yintercept = p_thr, color="black", linetype="dashed") +
  geom_hline(yintercept = -log10(0.05), color="black", linetype="solid") +
  labs(title = sprintf("Auto-antibodies and Herpesviruses in %s", cfg$title),
       y = expression('-log'[10]*'(P-value)'), color="Effect Size", shape="Direction") +
  theme_bw(base_size = 11) +
  theme(legend.position="bottom", axis.title.x=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        plot.title=element_text(hjust=0.5))
ggsave(file.path(opt$outdir, sprintf("manhattan_%s.png", cfg$title)), p_manh, width=12, height=5, dpi=300)

# ---- Circos (chord) with prevalence labels ----
prev_h <- fread(opt$prev_human)  # gene_symbol, Prevalence
prev_v <- fread(opt$prev_virus)  # taxon_species, Prevalence
prev_gene <- setNames(prev_h$Prevalence, prev_h$gene_symbol)
prev_spec <- setNames(prev_v$Prevalence, prev_v$taxon_species)

rep_summary <- dat %>%
  group_by(taxon_species, gene_symbol) %>%
  summarize(association_count = n(), min_p = min(ifelse(is.na(P.adj), P, P.adj), na.rm=TRUE), .groups="drop") %>%
  mutate(log10p = -log10(min_p))
rep_summary$taxon_species <- factor(rep_summary$taxon_species,
                                    levels=c("HSV-1","HSV-2","VZV","CMV","HHV-6A","HHV-6B","HHV-7","EBV","HHV-8"))

auto_list  <- unique(rep_summary$gene_symbol)
virus_list <- unique(as.character(rep_summary$taxon_species))
all_sectors <- c(auto_list, virus_list)

min_val <- round(min(rep_summary$log10p, na.rm=TRUE), 2)
max_val <- round(max(rep_summary$log10p, na.rm=TRUE), 2)
col_fun <- circlize::colorRamp2(seq(min_val, max_val, length.out = 9),
                                colorRampPalette(c("#FEE5D9", "#A50F15"))(9))

circlize::circos.clear()
circlize::circos.par(start.degree = 110,
                     gap.after = c(rep(3, length(auto_list)-1), 12, rep(3, length(virus_list)-1), 12),
                     circle.margin = 0.3, track.margin = c(0.05,0.05))

png(file.path(opt$outdir, sprintf("circos_%s.png", cfg$title)), width=1800, height=1800, res=250)
chord_data <- rep_summary %>% rename(from = gene_symbol, to = taxon_species)
circlize::chordDiagram(
  x = chord_data[, c("from","to","association_count")],
  order = all_sectors,
  col = col_fun(chord_data$log10p),
  transparency = 0.3,
  directional = 0,
  annotationTrack = "grid",
  link.lwd = pmax(1, chord_data$log10p / max_val * 3)
)
circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
  sector.index <- circlize::get.cell.meta.data("sector.index")
  circlize::circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 1.1, sector.index,
                        facing='clockwise', niceFacing=TRUE, adj=c(0,0.5), cex=0.7, font=2)
  prev_val <- ifelse(sector.index %in% names(prev_gene), prev_gene[[sector.index]], prev_spec[[sector.index]])
  if (!is.na(prev_val)) {
    circlize::circos.text(CELL_META$xcenter, CELL_META$ylim[1] - 0.1,
                          labels = sprintf("%.2f", prev_val),
                          facing="reverse.clockwise", niceFacing=TRUE, adj=c(0,0.5), cex=0.6, font=3)
  }
})
grid::grid.draw(ComplexHeatmap::Legend(title=expression('-log'[10]*'(P)'), col_fun=col_fun, at=seq(min_val,max_val,length.out=6)))
dev.off()


Rscript scripts/plots_qtl_panels.R \
  --cohort MGBB-LLF \
  --annot-rds results/annot/MGBB-LLF_glm_annotated.rds \
  --prev-human results/annot/MGBB-LLF_prevalence_human.csv \
  --prev-virus results/annot/MGBB-LLF_prevalence_virus_species.csv

Rscript scripts/plots_qtl_panels.R \
  --cohort MGBB-ABC \
  --annot-rds results/annot/MGBB-ABC_glm_annotated.rds \
  --prev-human results/annot/MGBB-ABC_prevalence_human.csv \
  --prev-virus results/annot/MGBB-ABC_prevalence_virus_species.csv


