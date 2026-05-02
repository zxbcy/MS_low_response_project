# 02_prepare_phenotype_GSE138064.R
# Clean phenotype metadata and detect complete/partial responder labels.
# Stops with a clear message if labels are ambiguous.

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(janitor)
})

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

eset_path <- "data/raw/GSE138064_eset.rds"
if (!file.exists(eset_path)) {
  stop("Missing data/raw/GSE138064_eset.rds. Run scripts/01_download_GSE138064.R first.")
}

eset <- readRDS(eset_path)
pheno <- Biobase::pData(eset) %>%
  tibble::rownames_to_column("sample_id") %>%
  clean_names()

# Collapse row text to help pattern scanning across all columns.
collapsed <- apply(pheno, 1, function(x) paste(x, collapse = " | "))
collapsed_low <- tolower(collapsed)

# Heuristic detection for responder categories.
is_complete <- str_detect(collapsed_low, "complete responder|complete response|responder: complete|response: complete|ifn-beta complete")
is_partial <- str_detect(collapsed_low, "partial responder|partial response|responder: partial|response: partial|ifn-beta partial")

group <- rep(NA_character_, nrow(pheno))
group[is_complete & !is_partial] <- "complete_responder"
group[is_partial & !is_complete] <- "partial_responder"

pheno_clean <- pheno %>%
  mutate(response_group = group,
         response_group_binary = case_when(
           response_group == "partial_responder" ~ 1L,
           response_group == "complete_responder" ~ 0L,
           TRUE ~ NA_integer_
         ))

write.csv(pheno_clean, "data/processed/GSE138064_pheno_clean_candidate.csv", row.names = FALSE)

n_complete <- sum(pheno_clean$response_group == "complete_responder", na.rm = TRUE)
n_partial <- sum(pheno_clean$response_group == "partial_responder", na.rm = TRUE)
n_missing <- sum(is.na(pheno_clean$response_group))

if (n_complete == 0 || n_partial == 0 || n_missing > 0) {
  stop(
    paste0(
      "Ambiguous or incomplete responder labels detected.\n",
      "Candidate metadata exported: data/processed/GSE138064_pheno_clean_candidate.csv\n",
      "Please manually confirm complete_responder vs partial_responder labels, then save:\n",
      "data/processed/GSE138064_pheno_final.csv\n",
      "Counts -> complete: ", n_complete,
      ", partial: ", n_partial,
      ", unlabeled: ", n_missing
    )
  )
}

write.csv(pheno_clean, "data/processed/GSE138064_pheno_final.csv", row.names = FALSE)
message("Phenotype labels detected successfully and saved to data/processed/GSE138064_pheno_final.csv")
