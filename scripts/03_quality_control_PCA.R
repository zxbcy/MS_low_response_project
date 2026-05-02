# 03_quality_control_PCA.R
# Perform expression preprocessing, optional log2 transform, gene mapping, deduplication, and PCA.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/PCA", recursive = TRUE, showWarnings = FALSE)

eset <- readRDS("data/raw/GSE138064_eset.rds")
pheno <- read.csv("data/processed/GSE138064_pheno_final.csv", stringsAsFactors = FALSE)

expr <- Biobase::exprs(eset)
feature <- Biobase::fData(eset) %>% tibble::rownames_to_column("probe_id")

# Log2 transform if needed based on quantile heuristic.
qx <- quantile(expr, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE)
need_log2 <- (qx[5] > 100) || ((qx[6] - qx[1]) > 50)
if (need_log2) {
  expr[expr <= 0] <- NA
  expr <- log2(expr)
}

expr_df <- as.data.frame(expr) %>% tibble::rownames_to_column("probe_id")

# Attempt to identify gene symbol column.
symbol_col_candidates <- c("gene_symbol", "gene.symbol", "symbol", "genesymbol", "gene assignment")
symbol_col <- names(feature)[tolower(names(feature)) %in% symbol_col_candidates]

if (length(symbol_col) == 0) {
  # Fallback: detect any column containing "symbol" keyword.
  symbol_col <- names(feature)[str_detect(tolower(names(feature)), "symbol")]
}

if (length(symbol_col) == 0) {
  stop("No gene symbol annotation column found in feature metadata.")
}

symbol_col <- symbol_col[1]

merged <- expr_df %>%
  left_join(feature %>% select(probe_id, gene_symbol = all_of(symbol_col)), by = "probe_id") %>%
  filter(!is.na(gene_symbol), gene_symbol != "")

mat <- as.matrix(merged[, setdiff(colnames(merged), c("probe_id", "gene_symbol"))])
rownames(mat) <- merged$probe_id

# Remove duplicate gene symbols by keeping probe with highest variance.
probe_var <- apply(mat, 1, var, na.rm = TRUE)
merged$probe_var <- probe_var

dedup <- merged %>%
  group_by(gene_symbol) %>%
  slice_max(order_by = probe_var, n = 1, with_ties = FALSE) %>%
  ungroup()

expr_gene <- as.matrix(dedup[, setdiff(colnames(dedup), c("probe_id", "gene_symbol", "probe_var"))])
rownames(expr_gene) <- dedup$gene_symbol

write.csv(expr_gene, "data/processed/GSE138064_expr_gene_dedup.csv", row.names = TRUE)

# PCA on samples.
expr_pca <- t(expr_gene)
pca <- prcomp(expr_pca, center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca$x[, 1:2]) %>% tibble::rownames_to_column("sample_id")

plot_df <- pca_df %>% left_join(pheno %>% select(sample_id, response_group), by = "sample_id")

p <- ggplot(plot_df, aes(PC1, PC2, color = response_group, label = sample_id)) +
  geom_point(size = 3, alpha = 0.9) +
  theme_bw() +
  labs(title = "PCA: GSE138064", subtitle = "Complete vs Partial responders", color = "Group")

ggsave("figures/PCA/GSE138064_PCA.png", p, width = 8, height = 6, dpi = 300)
dir.create("results/QC", recursive = TRUE, showWarnings = FALSE)
write.csv(plot_df, "results/QC/GSE138064_PCA_scores.csv", row.names = FALSE)
