# Project Plan

## Title
Risk stratification of low response to conventional therapy in multiple sclerosis by integrating peripheral blood transcriptomics, cerebrospinal fluid immune signatures, and MRI lesion burden.

## Initial implemented module
- Public transcriptomics dataset: GSE138064
- Comparison: IFN-beta complete responders vs partial responders

## Objectives
1. Build a reproducible preprocessing and QC pipeline.
2. Identify differential transcriptomic features for low-response biology.
3. Derive a compact predictive signature using LASSO logistic regression.
4. Prepare structure for future CSF and MRI integration.

## Milestones
1. Data acquisition and phenotype curation.
2. QC and exploratory analysis.
3. Differential expression and visualization.
4. Functional enrichment (GO/KEGG).
5. Classification modeling and ROC/AUC evaluation.
6. Extension planning for multi-modal integration.

## Reproducibility principles
- Relative paths only.
- Script-level directory creation.
- Explicit stopping on ambiguous phenotype labels.
- No large data/results artifacts committed to Git.
