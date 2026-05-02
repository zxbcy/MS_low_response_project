# 05_volcano_heatmap.R
# Generate volcano plot and heatmap of top differential genes.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

dir.create("figures/volcano", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/heatmap", recursive = TRUE, showWarnings = FALSE)

deg <- read.csv("results/DEG/GSE138064_limma_all_genes.csv", stringsAsFactors = FALSE)
expr_gene <- read.csv("data/processed/GSE138064_expr_gene_dedup.csv", row.names = 1, check.names = FALSE)
pheno <- read.csv("data/processed/GSE138064_pheno_final.csv", stringsAsFactors = FALSE)

deg$significance <- ifelse(deg$adj.P.Val < 0.05 & abs(deg$logFC) >= 1, "Significant", "Not Significant")

p_volcano <- ggplot(deg, aes(x = logFC, y = -log10(P.Value), color = significance)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey70")) +
  theme_bw() +
  labs(title = "Volcano plot", x = "log2 Fold Change", y = "-log10(P-value)")

top_labels <- head(deg[order(deg$adj.P.Val), ], 15)
p_volcano <- p_volcano +
  geom_text_repel(data = top_labels, aes(label = gene_symbol), size = 3, max.overlaps = 100)

ggsave("figures/volcano/GSE138064_volcano.png", p_volcano, width = 8, height = 6, dpi = 300)

# Heatmap of top 50 genes by adjusted p-value.
top_genes <- head(deg$gene_symbol[order(deg$adj.P.Val)], 50)
mat <- as.matrix(expr_gene[top_genes, , drop = FALSE])

annotation_col <- data.frame(response_group = pheno$response_group)
rownames(annotation_col) <- pheno$sample_id
annotation_col <- annotation_col[colnames(mat), , drop = FALSE]

pheatmap(mat,
         scale = "row",
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = FALSE,
         filename = "figures/heatmap/GSE138064_top50_heatmap.png",
         width = 9,
         height = 11)
