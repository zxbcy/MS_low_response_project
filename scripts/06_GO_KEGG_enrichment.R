# 06_GO_KEGG_enrichment.R
# Perform GO and KEGG enrichment with clusterProfiler.

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)

deg_sig <- read.csv("results/DEG/GSE138064_limma_sig_genes.csv", stringsAsFactors = FALSE)

if (nrow(deg_sig) == 0) {
  stop("No significant genes for enrichment. Check DEG thresholds.")
}

genes <- unique(deg_sig$gene_symbol)

gene_df <- bitr(genes,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

if (nrow(gene_df) == 0) {
  stop("No genes mapped to ENTREZID for enrichment.")
}

ego <- enrichGO(gene = gene_df$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                readable = TRUE)

kegg <- enrichKEGG(gene = gene_df$ENTREZID,
                   organism = "hsa",
                   pAdjustMethod = "BH")

write.csv(as.data.frame(ego), "results/enrichment/GSE138064_GO_BP.csv", row.names = FALSE)
write.csv(as.data.frame(kegg), "results/enrichment/GSE138064_KEGG.csv", row.names = FALSE)
