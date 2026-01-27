# MIPSA-HSV

HSV peptides and Human Antibodies in MIPSA
Analysis for: **Large-scale antibody reactome profiling identifies herpesvirus–autoantigen associations underlying chronic diseases**.

This repository reproduces the analyses described in the manuscript: cohort seroprevalence, Hu↔HSV GLMs, replication, prediction models, disease associations, and tripartite network plots.

## Data expectations
- **Viral**: VirSIGHT peptide/protein hits, with UniProt IDs and taxonomic annotations.
- **Human**: HuSIGHT full-length hits (binary), with UniProt and HGNC gene symbol.
- **Phenotypes**: Demographic covariates (Age, Sex/Gender, BMI, Race, Smoking, Alcohol, ICS, and CCI modified by age), and disease endpoints (prevalent, incident).
- **Replication cohort**: matching ID space for strict peptide–autoantigen ID matches.

## Repro order (discovery → replication → prediction → disease association)
### Statistics for diseases in the cohorts (Figure 1)
Rscript scripts/disease_counts.R \
  --llf_phe Data/llf_1290_phe.tsv \
  --lec_phe Data/lec_763_phe.tsv \
  --out_dir results/Figure1/

### Statistics for VirSIGHT and HuSIGHT (Figure 2)
Rscript scripts/seroprev_plots.R \
  --llf_promax Data/IB1007_VirSIGHT_Promax_Hits_Fold-Over-Background.csv \
  --abc_promax Data/IB1021_VirSIGHT_Promax_Hits_Fold-Over-Background.csv \
  --leo_promax Data/IB1189_VirSIGHT_Promax_Hits_Fold-Over-Background.csv \
  --llf_fl Data/IB1007_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv \
  --abc_fl Data/IB1021_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv \
  --leo_fl Data/IB1189_HuSIGHT_FullLength_Hits_Fold-Over-Background.tsv \
  --out_dir results/Figure2
   
### Build binary matrices from raw VirSIGHT/HuSIGHT CSVs
1) LLF cohort
Rscript scripts/build_binary_matrices.R \
  --cohort MGBB-LLF \
  --virsight_promax.llf Data/IB1007_VirSIGHT_Promax_Hits_Fold-Over-Background.csv \
  --husight_fl.llf Data/IB1007_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv \
  --pheno.llf Data/llf_1290_phe.tsv \
  --out_dir Data \

2) LEC cohort
Rscript scripts/build_binary_matrices.R \
  --cohort MGBB-LEC \
  --virsight_promax.abc Data/IB1021_VirSIGHT_Promax_Hits_Fold-Over-Background.csv \
  --husight_fl.abc Data/IB1021_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv \
  --virsight_promax.leo Data/IB1189_VirSIGHT_Promax_Hits_Fold-Over-Background_Old-Format-HSV.tsv \
  --husight_fl.leo Data/IB1189_HuSIGHT_FullLength_Hits_Fold-Over-Background.tsv \
  --pheno.lec Data/lec_763_phe.tsv \
  --out_dir Data \

### Hu–Virus pairwise GLMs (Hu Ab ~ Virus peptide + covariates)
Rscript scripts/pairwise_glm.R \
1) LLF cohort
Rscript scripts/run_hu_vs_hsv_glm.R \
  --llf_phe Data/llf_1290_phe.tsv \
  --llf_vi_promax Data/IB1007_VirSIGHT_Promax_Hits_Fold-Over-Background.csv \
  --llf_hu_fl Data/IB1007_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv \
  --llf_vi_bin Data/hsv_promax_bin_MGBB-LLF.tsv \
  --llf_hu_bin Data/human_fl_bin_MGBB-LLF.tsv \
  --out_dir results/Run \
  --n_cores 6

2) LEC cohort
  --llf_vi_promax Data/IB1007_VirSIGHT_Promax_Hits_Fold-Over-Background.csv \
  --llf_hu_fl Data/IB1007_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv \
  --lec_phe Data/lec_763_phe.tsv \
  --lec_vi_bin Data/hsv_promax_bin_MGBB-LEC.tsv \
  --lec_hu_bin Data/human_fl_bin_MGBB-LEC.tsv \
  --out_dir results/Run \
  --n_cores 6

### Manhattan, Diamond, and Circos plots (Figure 3)
Rscript scripts/plots_pairwise_panels.R \
  --llf_all results/Run/llf_hsv_bin_glm_all.tsv \
  --llf_sig results/Run/llf_hsv_bin_glm_sig.tsv \
  --llf_vi_promax Data/IB1007_VirSIGHT_Promax_Hits_Fold-Over-Background.csv \
  --llf_hu_fl Data/IB1007_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv \
  --llf_vi_bin Data/hsv_promax_bin_MGBB-LLF.tsv \
  --llf_hu_bin Data/human_fl_bin_MGBB-LLF.tsv \
  \
  --lec_test results/Run/lec_hsv_bin_glm_test.tsv \
  --lec_vi_bin Data/hsv_promax_bin_MGBB-LEC.tsv \
  --lec_hu_bin Data/human_fl_bin_MGBB-LEC.tsv \
  --rep_llf results/Run/hsv_bin_fchange_rep_llf.tsv \
  --rep_lec results/Run/hsv_bin_fchange_rep_lec.tsv \
  --out_dir results/Figure3
  
### Prediction (Figure 4): LASSO + ROC/PR + sens/spec + corr
Rscript scripts/lasso_prediction.R \
  --train_vi Data/hsv_promax_bin_MGBB-LLF.tsv \
  --train_hu Data/human_fl_bin_MGBB-LLF.tsv \
  --train_vi_anno Data/IB1007_VirSIGHT_Promax_Hits_Fold-Over-Background.csv \
  --train_hu_anno Data/IB1007_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv \
  --test_vi Data/hsv_promax_bin_MGBB-LEC.tsv \
  --test_hu Data/human_fl_bin_MGBB-LEC.tsv \
  --rep_hsv Run/hsv_bin_fchange_rep_llf.tsv \
  --out_dir results/Figure4/

### Disease GLMs (prevalent or incident disease)
Rscript scripts/Dx_prevalent_glm_hu_all.R
Rscript scripts/Dx_prevalent_glm_vi_all.R
Rscript scripts/Dx_incident_cox_hu_all.R
Rscript scripts/Dx_incident_cox_hu_all.R

  --llf_pheno Data/llf_1290_phe.tsv
  --llf_vi_promax Data/IB1007_VirSIGHT_Promax_Hits_Fold-Over-Background.csv
  --llf_hu_fl Data/IB1007_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv
  --llf_vi_bin Data/hsv_promax_bin_MGBB-LLF.tsv
  --llf_hu_bin Data/human_fl_bin_MGBB-LLF.tsv
  --rep_llf Run/hsv_bin_fchange_rep_llf.tsv
  --llf_dx_count results/Figure1/llf_case_counts_filtered.csv
  --lec_pheno Data/lec_763_phe.tsv
  --lec_vi_bin Data/hsv_promax_bin_MGBB-LEC.tsv
  --lec_hu_bin Data/human_fl_bin_MGBB-LEC.tsv
  --abc_vi_promax Data/IB1021_VirSIGHT_Promax_Hits_Fold-Over-Background.csv
  --abc_hu_fl Data/IB1021_HuSIGHT_FullLength_Hits_Fold-Over-Background.csv
  --leo_vi_promax Data/IB1189_VirSIGHT_Promax_Hits_Fold-Over-Background_Old-Format-HSV.tsv
  --leo_hu_fl Data/IB1189_HuSIGHT_FullLength_Hits_Fold-Over-Background.tsv
  --rep_lec Run/hsv_bin_fchange_rep_lec.tsv
  --lec_dx_count results/Figure1/lec_case_counts_filtered.csv
  --sig_ab results/Figure4/lec_test_auc_summary_sig.csv
  --out_dir results/Disease

### Network plots (Figure 5)
Rscript scripts/network_analysis.R
  --human_dx_inc results/Disease/rep_hu_dx_all_inc_coxphf_ann.tsv
  --virus_dx_inc results/Disease/rep_vi_dx_all_inc_coxphf_ann.tsv
  --human_dx_pre results/Disease/rep_hu_dx_all_pre_glmf_ann.tsv
  --virus_dx_pre results/Disease/rep_vi_dx_all_pre_glmf_ann.tsv
  --sig_ab Run/hsv_bin_fchange_rep_llf.tsv
  --out_dir results/Network
   
## Software
R ≥ 4.2 with: "argparse", "broom", "corrplot", "coxphf", "data.table", "doParallel", "dplyr", "foreach", "forcats", "ggraph", "ggplot2", "ggrepel", "glmnet", "grid", "hrbrthemes", "igraph", "logistf", "optparse", "parallel", "pROC", "purrr", "readr", "scales", "stringr", "survival", "tidyr", and "viridis".

## Citation
Please cite the manuscript when using this code.
