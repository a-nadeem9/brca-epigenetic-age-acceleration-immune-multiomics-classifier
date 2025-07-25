# BRCA Epigenetic Age Acceleration + Immune + Multi‑omics Classifier

> **Goal**: Quantify **epigenetic age acceleration** in breast cancer, integrate it with **immune scores** and other **multi‑omics** features, and train/evaluate machine‑learning classifiers (binary & multiclass). Epigenetic clocks (e.g., Horvath, PhenoAge, GrimAge) have been linked to breast cancer biology and risk, but results are mixed and tissue/clock‑dependent, motivating a careful, reproducible analysis like this one. [1–6]

<p align="left">
  <img alt="Python" src="https://img.shields.io/badge/Python-3.10+-blue.svg">
  <img alt="R" src="https://img.shields.io/badge/R-%3E=4.0-276DC3.svg">
  <img alt="Jupyter" src="https://img.shields.io/badge/Jupyter-%F0%9F%93%9A-orange.svg">
</p>

## Quick start

```bash
# 1) Clone
git clone https://github.com/a-nadeem9/brca-epigenetic-age-acceleration-immune-multiomics-classifier.git
cd brca-epigenetic-age-acceleration-immune-multiomics-classifier

# 2) Create & activate an environment
conda create -n brca-eaa python=3.10 -y
conda activate brca-eaa

# 3) Install Python deps
pip install -r requirements.txt

# 4) Launch notebooks
jupyter lab
```

> **R side:** The `scripts/` folder contains lightweight R utilities. Make sure you have **R ≥ 4.0** and install the required packages.



## Repository structure

```
brca-epigenetic-age-acceleration-immune-multiomics-classifier
├── notebooks
│   ├── 01a_data_preprocessing.ipynb
│   ├── 01b_rna_preprocessing.ipynb
│   ├── 01c_methylation_preprocessing.ipynb
│   ├── 02_delta_age_calculation.ipynb
│   ├── 03_immune_score.ipynb
│   ├── 04_binary_classification.ipynb
│   ├── 05_multiclass_classification.ipynb
│   └── 06_annotation.ipynb
├── scripts
│   ├── 01_download_tcga.R
│   ├── 02_merge_matched_data.R
│   ├── 03_tcga_subtype.R
│   └── 04_epigenetic_age.R
├── requirements.txt
└── README.md

```



## Notebook guide

| Notebook | What it does |
|---|---|
| **01a_data_preprocessing.ipynb** | Load/clean clinical & sample metadata, harmonise IDs, basic QC. |
| **01b_rna_preprocessing.ipynb** | RNA‑seq normalisation, filtering, optional batch correction; produce expression features. |
| **01c_methylation_preprocessing.ipynb** | Probe/array QC, normalisation, batch correction of methylation data. |
| **02_delta_age_calculation.ipynb** | Compute DNAm age via chosen clock(s) (e.g., Horvath, PhenoAge, GrimAge) and derive **ΔAge**. |
| **03_immune_score.ipynb** | Generate/ingest immune deconvolution scores (e.g., CIBERSORT, TIMER, xCell) and assemble immune features. |
| **04_binary_classification.ipynb** | Train & evaluate binary models (e.g., Basal vs. Luminal): AUROC/AUPRC, SHAP/feature importances. |
| **05_multiclass_classification.ipynb** | Train & evaluate multiclass models (e.g., PAM50): macro/micro metrics, confusion matrices, SHAP. |
| **06_annotation.ipynb** | Biological/clinical interpretation, pathway enrichment, figure/table export. |


## References

1. Karim et al. **Epigenetic Ageing and Breast Cancer Risk: A Systematic Review.** *Cancer Medicine* (2024).  
2. Horvath S. **DNA methylation age of human tissues and cell types.** *Genome Biology* (2013).  
3. Levine ME et al. **An epigenetic biomarker of aging for lifespan and healthspan (DNAm PhenoAge).** *Aging (Albany NY)* (2018).  
4. Lu AT et al. **DNA methylation GrimAge strongly predicts lifespan and healthspan.** *Aging (Albany NY)* (2019).  
5. Starnawska A et al. **Evaluation of epigenetic age acceleration scores and their associations with phenotypes.** *Genes* (2023).  
6. Snell LM et al. **Evidence of accelerated epigenetic aging of breast tissues in breast cancer.** *Breast Cancer Research* (2022).
