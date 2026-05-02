# 07_LASSO_signature.R
# Build LASSO logistic classifier and evaluate ROC/AUC.

suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
  library(ggplot2)
})

dir.create("results/models", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/ROC", recursive = TRUE, showWarnings = FALSE)

expr_gene <- read.csv("data/processed/GSE138064_expr_gene_dedup.csv", row.names = 1, check.names = FALSE)
pheno <- read.csv("data/processed/GSE138064_pheno_final.csv", stringsAsFactors = FALSE)

common_samples <- intersect(colnames(expr_gene), pheno$sample_id)
X <- t(as.matrix(expr_gene[, common_samples, drop = FALSE]))
ph <- pheno[match(common_samples, pheno$sample_id), ]
y <- ph$response_group_binary

if (any(is.na(y))) stop("response_group_binary contains NA values.")

set.seed(123)
cvfit <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 5, type.measure = "auc")
fit <- glmnet(X, y, family = "binomial", alpha = 1, lambda = cvfit$lambda.min)

pred <- as.numeric(predict(cvfit, newx = X, s = "lambda.min", type = "response"))
roc_obj <- roc(response = y, predictor = pred, quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

coef_df <- as.matrix(coef(fit))
write.csv(coef_df, "results/models/GSE138064_lasso_coefficients.csv")
write.csv(data.frame(sample_id = common_samples, y = y, pred = pred),
          "results/models/GSE138064_lasso_predictions.csv",
          row.names = FALSE)

roc_df <- data.frame(
  specificity = roc_obj$specificities,
  sensitivity = roc_obj$sensitivities
)

p <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(size = 1.1, color = "#2C7BB6") +
  geom_abline(linetype = "dashed", color = "grey50") +
  theme_bw() +
  labs(title = paste0("ROC Curve (AUC = ", round(auc_val, 3), ")"),
       x = "False Positive Rate",
       y = "True Positive Rate")

ggsave("figures/ROC/GSE138064_LASSO_ROC.png", p, width = 7, height = 6, dpi = 300)
