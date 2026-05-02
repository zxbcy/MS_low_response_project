# 04_differential_expression_limma.R
# Differential expression with limma for partial vs complete responders.

suppressPackageStartupMessages({
  library(limma)
})

dir.create("results/DEG", recursive = TRUE, showWarnings = FALSE)

expr_gene <- read.csv("data/processed/GSE138064_expr_gene_dedup.csv", row.names = 1, check.names = FALSE)
pheno <- read.csv("data/processed/GSE138064_pheno_final.csv", stringsAsFactors = FALSE)

# Ensure sample order consistency.
common_samples <- intersect(colnames(expr_gene), pheno$sample_id)
expr_gene <- as.matrix(expr_gene[, common_samples, drop = FALSE])
pheno <- pheno[match(common_samples, pheno$sample_id), ]

if (any(is.na(pheno$response_group))) stop("Missing response_group in phenotype table.")

group <- factor(pheno$response_group, levels = c("complete_responder", "partial_responder"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr_gene, design)
contrast <- makeContrasts(partial_responder - complete_responder, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

deg <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "P")
deg$gene_symbol <- rownames(deg)
write.csv(deg, "results/DEG/GSE138064_limma_all_genes.csv", row.names = FALSE)

sig <- subset(deg, adj.P.Val < 0.05 & abs(logFC) >= 1)
write.csv(sig, "results/DEG/GSE138064_limma_sig_genes.csv", row.names = FALSE)

message("DEG analysis complete.")
