# Multivariate GWAS Statistical Analysis Pipeline
**CDC 1.0.0** — Post-REGENIE multivariate module 

Author: Nadeem Khan · INRS-Centre Armand-Frappier Santé-Biotechnologie 
GitHub: [nkhan119](https://github.com/nkhan119)

---

## Overview

This pipeline utilized the output of `gwas-pipeline-v3` and runs four
analysis phases on all harmonised GWAS summary statistics
(`*_sumstats.tsv.gz` from Cohort_A and Cohort_B):

| Phase | Method | Key outputs |
|-------|--------|-------------|
| 1 · Hypothesis Testing | Pearson/Spearman cross-trait correlations, Hotelling T², Fisher Z CI | Correlation heatmaps, cross-cohort effect-size tornado plots |
| 2 · Regression | OLS, Ridge, Lasso, ElasticNet (α via 5-fold CV), PCR, PLS, VIF | R² comparison chart, VIF heatmap, model table |
| 3 · Multiple Testing | Bonferroni, Holm, BH-FDR, BY-FDR, Storey q-value, λ_GC | QQ-plots, λ_GC bubble chart, hit-count comparison |
| 4 · Causal Inference | 2SLS/IVW, MR-Egger intercept, Cochran Q/I², LOO sensitivity, causal network | Forest plots, network heatmap, LOO boxplot |

---

## Directory Structure

```
multivariate-pipeline/
├── main.nf                   
├── nextflow.config
├── modules/
│   └── multivariate_analysis.nf
├── bin/
│   └── multivariate_analysis.py   
├── envs/
│   └── multivariate_env.yaml      
├── docker/
│   └── Dockerfile                 
└── results/                        
    └── multivariate/
        ├── Multivariate_Report.html
        ├── phase1_hypothesis/
        ├── phase2_regression/
        ├── phase3_mtc/
        └── phase4_causal/
```

---

## Quick Start

### Option A

```bash
# 1. Create the environment once
conda env create -f envs/multivariate_env.yaml
conda activate cdc_multivariate

# 2. Run
nextflow run main.nf \
    -profile local,conda \
    --sumstats_dir ~/Downloads/1000G_phase3_common_norel/GWAS_analysis \
    --out_dir       results/ \
    -resume
```

### Option B

```bash
# 1. Build the image FROM THE REPO ROOT (not from docker/ subdirectory)
cd ~/Downloads/multivariate-pipeline
docker build -f docker/Dockerfile -t nkhan119/cdc-multivariate:1.0.0 .

# 2. Run
nextflow run main.nf \
    -profile docker \
    --sumstats_dir ~/Downloads/1000G_phase3_common_norel/GWAS_analysis \
    --out_dir       results/ \
    -resume
```

### Option C

```bash
# Build and push image to Docker Hub first
docker build -f docker/Dockerfile -t nkhan119/cdc-multivariate:1.0.0 .
docker push nkhan119/cdc-multivariate:1.0.0

# Then on slurm:
nextflow run main.nf \
    -profile narval,singularity \
    --sumstats_dir ~/projects/def-fveyrier/nad119/GWAS/1000G/GWAS_analysis \
    --out_dir       ~/projects/def-fveyrier/nad119/multivariate_results \
    -resume
```

---

## Integration with gwas-pipeline-v3

The `--sumstats_dir` should point to the **same directory** used as
`--out_dir` (or `--raw_dir`) in `gwas-pipeline-v3`. The pipeline
auto-discovers all files matching:

```
{sumstats_dir}/**/sumstats/*_sumstats.tsv.gz
```

This covers the gwas-pipeline-v3 layout:
```
GWAS_analysis/
  gwas/Cohort_A/sumstats/ldl_cholesterol_sumstats.tsv.gz
  gwas/Cohort_B/sumstats/bmi_sumstats.tsv.gz
  ...
```

---

## Expected Sumstat Columns

Minimum required:

| Column | Description |
|--------|-------------|
| `SNP`  | Variant ID (chr:pos or rsID) |
| `CHR`  | Chromosome |
| `BP`   | Base-pair position |
| `BETA` | Effect size |
| `SE`   | Standard error |
| `P`    | p-value |

Optional (used if present): `Z`, `A1FREQ` (MAF), `INFO`, `A1`, `A2`

---

## Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--sumstats_dir` | (required) | Directory tree containing `*_sumstats.tsv.gz` |
| `--out_dir` | `./results` | Output directory |
| `--author` | `Nadeem Khan` | Report author name |
| `--institute` | INRS-CAFSB | Report institute |
| `--docker_image` | `nkhan119/cdc-multivariate:1.0.0` | Docker image tag |

---

## Outputs

```
results/multivariate/
├── Multivariate_Report.html
│
├── phase1_hypothesis/
│   ├── cross_trait_correlations.tsv      
│   ├── hotelling_T2.tsv                  
│   ├── Phase1_HypothesisTesting.xlsx     
│   └── static/
│       ├── corr_heatmap_Cohort_A.pdf/.png
│       ├── corr_heatmap_Cohort_B.pdf/.png
│       ├── tornado_ldl_cholesterol.pdf/.png
│       ├── tornado_bmi.pdf/.png           
│       └── hotelling_T2.pdf/.png
│
├── phase2_regression/
│   ├── regression_comparison.tsv        
│   ├── vif_diagnostics.tsv              
│   ├── Phase2_Regression.xlsx        
│   └── static/
│       ├── r2_comparison.pdf/.png
│       ├── vif_heatmap.pdf/.png
│       └── regularisation_r2.pdf/.png
│
├── phase3_mtc/
│   ├── multiple_testing_correction.tsv
│   ├── genomic_inflation.tsv
│   ├── Phase3_MultiTesting.xlsx
│   └── static/
│       ├── qq_Cohort_A_ldl_cholesterol.pdf/.png
│       ├── qq_Cohort_B_bmi.pdf/.png
│       ├── lambda_gc_bar.pdf/.png
│       ├── sig_hits_grouped.pdf/.png
│       └── storey_pi0.pdf/.png
│
└── phase4_causal/
    ├── iv_causal_estimates.tsv
    ├── loo_sensitivity.tsv
    ├── Phase4_CausalInference.xlsx
    └── static/
        ├── forest_Cohort_A.pdf/.png
        ├── forest_Cohort_B.pdf/.png
        ├── causal_network_Cohort_A.pdf/.png
        ├── causal_network_Cohort_B.pdf/.png
        ├── loo_Cohort_A.pdf/.png
        ├── loo_Cohort_B.pdf/.png
        ├── funnel_Cohort_A.pdf/.png
        └── funnel_Cohort_B.pdf/.png
```

---

## Building the Docker Image

```bash
cd docker/
docker build -t nkhan119/cdc-multivariate:1.0.0 .

# Optional: push to Docker Hub
docker push nkhan119/cdc-multivariate:1.0.0
```

The image is based on `python:3.11-slim-bullseye` and includes:
`numpy`, `scipy`, `pandas`, `statsmodels`, `scikit-learn`,
`plotly`, `kaleido`, `dowhy`, `econml`, `pingouin`, `lifelines`
