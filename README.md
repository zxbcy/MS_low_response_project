# Risk stratification of low response to conventional therapy in multiple sclerosis

This repository contains **code and documentation only** for a reproducible transcriptomics-first workflow focused on:

> Risk stratification of low response to conventional therapy in multiple sclerosis by integrating peripheral blood transcriptomics, cerebrospinal fluid immune signatures, and MRI lesion burden.

## Scope of this initial module

Current implementation covers the first transcriptomics module for **GSE138064** (peripheral blood), with primary comparison:

- **IFN-beta complete responders vs partial responders**

Biological question:

- Identify peripheral blood transcriptomic features associated with low response to conventional therapy in multiple sclerosis.

## Repository policy (important)

- Keep this repository to **code and documentation only**.
- Do **not** commit raw GEO files, MRI NIfTI files, processed matrices, result tables, model files, or generated figures.
- Use **relative paths only** (no hard-coded `C:\` or other absolute drive paths).

## Directory structure

- `scripts/`: analysis scripts
- `docs/`: planning and data-source documentation
- runtime-created folders (ignored by git):
  - `data/raw`
  - `data/processed`
  - `results/DEG`
  - `results/enrichment`
  - `results/models`
  - `figures/PCA`
  - `figures/volcano`
  - `figures/heatmap`
  - `figures/ROC`

## Prerequisites

Install R (>=4.2 recommended) and the required packages.

### 1) Install CRAN packages

```r
install.packages(c(
  "tidyverse", "pheatmap", "ggrepel", "glmnet", "pROC", "janitor"
))
```

### 2) Install Bioconductor packages

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "GEOquery", "limma", "Biobase", "clusterProfiler", "org.Hs.eg.db", "AnnotationDbi"
))
```

## Run order (step-by-step)

From the repository root (for example: `D:\MS_low_response_project`), run scripts in order:

```r
source("scripts/01_download_GSE138064.R")
source("scripts/02_prepare_phenotype_GSE138064.R")
source("scripts/03_quality_control_PCA.R")
source("scripts/04_differential_expression_limma.R")
source("scripts/05_volcano_heatmap.R")
source("scripts/06_GO_KEGG_enrichment.R")
source("scripts/07_LASSO_signature.R")
```

## Notes on responder labels

Script `02_prepare_phenotype_GSE138064.R` attempts to parse GEO metadata and detect responder labels.

- If complete responder / partial responder labels are confidently identified, it exports a clean phenotype table.
- If labels are ambiguous, it writes a candidate metadata table and **stops with a clear message** asking for manual confirmation.

This protects downstream analyses from silent label misclassification.

## Outputs

All generated outputs are written to ignored local folders via relative paths, including:

- cleaned phenotype metadata
- DEG results from limma
- PCA / volcano / heatmap / ROC figures
- GO and KEGG enrichment tables
- LASSO model artifacts

