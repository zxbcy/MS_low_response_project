# Data Sources

## Public source (current module)
- GEO accession: **GSE138064**
- Access path: downloaded programmatically via `GEOquery::getGEO()`

## Data handling policy
- Raw downloads are stored locally under `data/raw/`.
- Cleaned metadata and expression matrices are stored locally under `data/processed/`.
- No raw/processed data are committed to this repository.

## Future planned sources
- CSF immune signatures (cohort-specific or public references TBD)
- MRI lesion burden metrics (site-specific preprocessing pipeline TBD)

## Compliance notes
- Repository is restricted to code and documentation.
- All data-heavy outputs are excluded through `.gitignore`.
