if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages(c("tidyverse", "limma", "edgeR"))
library(tidyverse)
library(edgeR)
library(readr)
library(stringr)
library(purrr)

# === Define paths ===
base_path <- "D:/Research/GEO/GSE146853"
counts_file <- file.path(base_path, "raw_counts.csv")
metadata_file <- file.path(base_path, "metadata.csv")
filtered_counts_file <- file.path(base_path, "filtered_counts.csv")
filtered_metadata_file <- file.path(base_path, "filtered_metadata.csv")
results_path <- file.path(base_path, "Results", "samples_with_replications")
out_path <- file.path(results_path, "EdgeR")
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

# === Preprocess counts and metadata ===
# Metadata
metadata <- read.csv(metadata_file, row.names = "SampleName")
metadata$Sample_ID <- trimws(metadata$Sample_ID)
metadata <- metadata %>%
  mutate(Sample_ID = gsub("-", ".", gsub("_GeneCount", "", Sample_ID))) %>%
  filter(!is.na(Disease), Disease != "NA", Disease != "")
metadata$Patient_ID <- sub("\\.\\d+$", "", metadata$Sample_ID)

replicated <- metadata %>%
  group_by(Patient_ID) %>%
  tally() %>%
  filter(n == 2) %>%
  pull(Patient_ID)

filtered_metadata <- metadata %>%
  filter(Patient_ID %in% replicated)

# Counts
raw_counts <- read.csv(counts_file, row.names = "GeneID")
colnames(raw_counts) <- colnames(raw_counts) %>%
  trimws() %>%
  gsub("_GeneCount", "", .) %>%
  gsub("-", ".", .)

common_samples <- intersect(colnames(raw_counts), filtered_metadata$Sample_ID)
if (length(common_samples) == 0) {
  stop("No matching sample names found between raw_counts and metadata.")
}

filtered_counts <- raw_counts[, common_samples]
filtered_metadata <- filtered_metadata %>%
  filter(Sample_ID %in% common_samples) %>%
  arrange(match(Sample_ID, common_samples))

dge <- DGEList(counts = filtered_counts, group = filtered_metadata$Disease)
keep <- filterByExpr(dge)
write.csv(dge$counts[keep, ], filtered_counts_file)
write.csv(filtered_metadata, filtered_metadata_file)

# === edgeR Differential Expression ===
patients <- filtered_metadata %>%
  filter(Disease != "Healthy") %>%
  distinct(Patient_ID, Disease)

counts <- read.csv(filtered_counts_file, row.names = 1)
metadata <- read.csv(filtered_metadata_file, row.names = 1)

run_edger <- function(patient_id, disease, samples) {
  group <- factor(ifelse(samples$Disease == "Healthy", "Healthy", disease), levels = c("Healthy", disease))
  dge <- DGEList(counts = counts[, samples$Sample_ID], group = group) %>%
    calcNormFactors() %>%
    { .[filterByExpr(., group), , keep.lib.sizes = FALSE] }
  design <- model.matrix(~group)
  fit <- glmQLFit(estimateDisp(dge, design), design)
  res <- topTags(glmQLFTest(fit, coef = 2), n = Inf)$table %>%
    rownames_to_column("GeneID")
  write_csv(res, file.path(out_path, sprintf("edgeR_%s_%s.csv", patient_id, disease)))
}

for (i in seq_len(nrow(patients))) {
  patient_id <- patients$Patient_ID[i]
  disease <- patients$Disease[i]
  samples <- metadata %>% filter(Patient_ID == patient_id | Disease == "Healthy")
  run_edger(patient_id, disease, samples)
}

# === Merge all edgeR results (significant genes only) ===
edgeR_files <- list.files(out_path, full.names = TRUE, pattern = "^edgeR_s_.*\\.csv$")

merge_edger_results <- function(file_path) {
  df <- read_csv(file_path, show_col_types = FALSE)
  df <- df %>%
    filter(PValue < 0.05) %>%
    mutate(
      Patient = str_match(basename(file_path), "^edgeR_s_(.*?)_")[,2],
      Disease = str_match(basename(file_path), "^edgeR_s_.*?_(.*?)\\.csv$")[,2]
    )
  return(df)
}

merged_edgeR <- map_dfr(edgeR_files, merge_edger_results)
write_csv(merged_edgeR, file.path(out_path, "edgeR_filtered_merged.csv"))
