# MIPSA-HSV

HSV peptides and Human Antibodies in MIPSA
Analysis for: **Large-scale antibody reactome profiling identifies herpesvirus–autoantigen associations underlying chronic diseases**.

This repository reproduces the analyses described in the manuscript: cohort seroprevalence, Hu↔HSV GLMs, replication, prediction models, disease associations, and tripartite network plots.

## Data expectations
- **Viral**: VirSIGHT peptide/protein hits (binary), with UniProt IDs and taxonomic annotations.
- **Human**: HuSIGHT full-length hits (binary), with UniProt and HGNC gene symbol.
- **Phenotypes**: Demographic covariates (Age, Sex/Gender, BMI, Race, Smoking, Alcohol), and disease endpoints (prevalent, incident).
- **Replication cohort**: matching ID space for strict peptide–autoantigen ID matches.

## Repro order (discovery → replication → prediction → disease association)
1) Statistics for diseases in the cohorts (Figure 1)
Rscript scripts/disease_counts.R \
  --llf_phe /data/MIPSA_Asthma_1290.csv \
  --abc_phe /data/Ab_pheno.csv \
  --out_dir results/Figure1/

2) Statistics for VirSIGHT and HuSIGHT (Figure 2)
Rscript scripts/virus_seroprev_plots.R \
  --llf_promax /data/IB1007_VirSIGHT_Promax_Hits_Fold-Over-Background.csv \
  --llf_varscore /data/IB1007_VirSIGHT_VARscores.csv \
  --abc_promax /data/IB1021_VirSIGHT_Promax_Hits_Fold-Over-Background.csv \
  --abc_varscore /data/IB1021_VirSIGHT_VARscores.csv \
  --out_dir results/Figure2 
Rscript scripts/human_seroprev_ridgeline.R \
  --llf_hu_fl /data/IB1007_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv \
  --abc_hu_fl /data/IB1021_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv \
  --out_dir results/Figure2
   
3) Build binary matrices from raw VirSIGHT/HuSIGHT CSVs
LLF cohort
Rscript scripts/build_binary_matrices.R \
  --cohort MGBB-LLF \
  --virsight_promax /data/IB1007_VirSIGHT_Promax_Hits_Fold-Over-Background.csv \
  --husight_fl      /data/IB1007_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv \
  --out_dir results/binaries \
  --min_prev_pct 1

ABC cohort
Rscript scripts/build_binary_matrices.R \
  --cohort MGBB-ABC \
  --virsight_promax /data/IB1021_VirSIGHT_Promax_Hits_Fold-Over-Background.csv \
  --husight_fl      /data/IB1021_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv \
  --out_dir results/binaries \
  --min_prev_pct 1


4) Hu–Virus pairwise GLMs (Hu Ab ~ Virus peptide + covariates), chunkable
Rscript scripts/run_hu_vs_hsv_glm.R \
  --phe data/MIPSA_Asthma_1290.csv \
  --vi_bin results/virus_proteins_binary.txt \
  --hu_bin results/human_fl_binary.txt \
  --vi_anno data/IB1007_VirSIGHT_Promax_Hits_Fold-Over-Background.csv \
  --rep_sig data/rep_hsv.sig \
  --aa 1 --bb 2258 \
  --out_dir results/Res_hu_vs_vi/
   
5) Collect + annotate + FDR + replication flag
Rscript scripts/collect_annotate_results.R \
  --res_dir results/Res_hu_vs_vi/ \
  --vi_anno data/IB1007_VirSIGHT_Promax_Hits_Fold-Over-Background.csv \
  --hu_anno data/IB1007_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv \
  --rep_sig data/rep_hsv.sig \
  --out_prefix results/qtl_hu_vs_vi

6) Manhattan, Diamond, and Circos plots (Figure 3)
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

7) Disease GLMs (combined prevalent + incident)
Rscript scripts/run_disease_glm.R \
  --phe data/MIPSA_Asthma_1290.csv \
  --vi_bin results/hsv_proteins_binary.txt \
  --hu_bin results/human_fl_binary.txt \
  --dx_prev data/dx_prevalence_sums \
  --dx_inc data/dx_incidence_sums \
  --out_dir results/Res_disease/

8) Prediction (Figure 4): LASSO + ROC/PR + sens/spec + corr
Rscript scripts/figure4_lasso_prediction.R \
  --train_vi results/virus_proteins_binary.txt \
  --train_hu results/human_fl_binary.txt \
  --train_vi_anno data/IB1007_VirSIGHT_Promax_Hits_Fold-Over-Background.csv \
  --train_hu_anno data/IB1007_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv \
  --test_vi data/ABC_IB1021_VirSIGHT_binary.txt \
  --test_hu data/ABC_IB1021_HuSIGHT_binary.txt \
  --ab_list PHLDA1:h6808,ZNF550:h14835,IQCB1:h3143,DNAJC12:h8069,P3H4:h8250 \
  --species_map "Human alphaherpesvirus 1=HSV-1;Human betaherpesvirus 5=CMV;Human gammaherpesvirus 4=EBV" \
  --out_dir results/Figure4/

9) Network plots (Figure 5)
Rscript scripts/network_plots.R \
  --human_dx results/rep_bin_fchange_FL_dx_com_glm_ann.txt \
  --proid_dx results/rep_bin_fchange_HSV_dx_com_glm_ann.txt \
  --rep_sig data/rep_hsv.sig \
  --virus "Human alphaherpesvirus 1" \
  --out_png results/network_plot_HSV1.png
   
## Software
R ≥ 4.2 with: `data.table`, `dplyr`, `tidyr`, `ggplot2`, `ggrepel`, `glmnet`, `pROC`, `igraph`, `ggraph`, `scales`, `argparse`, `readr`, `stringr`.

## Citation
Please cite the manuscript when using this code.
