# 01_download_GSE138064.R
# Download GSE138064 from GEO using GEOquery and save expression + metadata locally.

suppressPackageStartupMessages({
  library(GEOquery)
})

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

message("Downloading GEO series GSE138064 ...")
gse_list <- getGEO("GSE138064", GSEMatrix = TRUE, getGPL = TRUE)

if (length(gse_list) < 1) {
  stop("No ExpressionSet returned for GSE138064.")
}

# Use the first platform by default; downstream scripts keep platform metadata.
eset <- gse_list[[1]]
saveRDS(eset, file = "data/raw/GSE138064_eset.rds")

pheno <- Biobase::pData(eset)
feature <- Biobase::fData(eset)
expr <- Biobase::exprs(eset)

write.csv(pheno, "data/processed/GSE138064_pheno_raw.csv", row.names = TRUE)
write.csv(feature, "data/processed/GSE138064_feature_raw.csv", row.names = TRUE)
write.csv(expr, "data/processed/GSE138064_expr_raw.csv", row.names = TRUE)

message("Saved raw GEO objects and extracted tables to data/raw and data/processed.")
