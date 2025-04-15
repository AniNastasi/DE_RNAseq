if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages(c("tidyverse", "limma", "edgeR", "DESeq2"))
library(tidyverse)
library(limma)
library(edgeR)
library(DESeq2)
library(readr)
library(stringr)
library(purrr)

# === Define paths ===
base_path <- "D:/Research/GEO/GSE146853"
counts_file <- file.path(base_path, "raw_counts.csv")
metadata_file <- file.path(base_path, "metadata.csv")
filtered_counts_file <- file.path(base_path, "filtered_counts.csv")
filtered_metadata_file <- file.path(base_path, "filtered_metadata.csv")
results_path <- file.path(base_path, "Results")
out_paths <- list(DESeq2 = "DESeq2", EdgeR = "EdgeR", Limma = "limma", Merged = "Merged")
walk(out_paths, ~ dir.create(file.path(results_path, .x), showWarnings = FALSE))

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
head(raw_counts)
colnames(raw_counts) <- colnames(raw_counts) %>%
  trimws() %>%
  gsub("_GeneCount", "", .) %>%
  gsub("-", ".", .)

common_samples <- intersect(colnames(raw_counts), filtered_metadata$Sample_ID)
if (length(common_samples) == 0) {
  stop("No matching sample names found between raw_counts and metadata.")
} # Check if any common samples exist

filtered_counts <- raw_counts[, common_samples]
filtered_metadata <- filtered_metadata %>%
  filter(Sample_ID %in% common_samples) %>%
  arrange(match(Sample_ID, common_samples))

dge <- DGEList(counts = filtered_counts, group = filtered_metadata$Disease)
keep <- filterByExpr(dge)
write.csv(dge$counts[keep, ], filtered_counts_file)
write.csv(filtered_metadata, filtered_metadata_file)


# === Differential Expression Functions ===
patients <- filtered_metadata %>%
  filter(Disease != "Healthy") %>%
  distinct(Patient_ID, Disease)

counts <- read.csv(filtered_counts_file, row.names = 1)
metadata <- read.csv(filtered_metadata_file, row.names = 1)

# DESeq2
run_deseq2 <- function(patient_id, disease, samples) {
  cond <- ifelse(samples$Disease == "Healthy", "Healthy", disease)
  dds <- DESeqDataSetFromMatrix(round(counts[, samples$Sample_ID]), 
                                data.frame(row.names = samples$Sample_ID, Condition = factor(cond, levels = c("Healthy", disease))), 
                                design = ~Condition)
  res <- results(DESeq(dds)) %>% as.data.frame()
  write.csv(res, file.path(results_path, "DESeq2", sprintf("DESeq_%s_%s.csv", patient_id, disease)))
}

# edgeR
run_edger <- function(patient_id, disease, samples) {
  group <- factor(ifelse(samples$Disease == "Healthy", "Healthy", disease), levels = c("Healthy", disease))
  dge <- DGEList(counts = counts[, samples$Sample_ID], group = group) %>%
    calcNormFactors() %>%
    { .[filterByExpr(., group), , keep.lib.sizes = FALSE] }
  design <- model.matrix(~group)
  fit <- glmQLFit(estimateDisp(dge, design), design)
  res <- topTags(glmQLFTest(fit, coef = 2), n = Inf)$table %>%
    rownames_to_column("GeneID")
  write_csv(res, file.path(results_path, "EdgeR", sprintf("edgeR_%s_%s.csv", patient_id, disease)))
}

# limma+voom
run_limma <- function(patient_id, disease, samples) {
  group <- factor(ifelse(samples$Disease == "Healthy", "Healthy", disease), levels = c("Healthy", disease))
  dge <- DGEList(counts = counts[, samples$Sample_ID]) %>%
    calcNormFactors() %>%
    { .[filterByExpr(., group), , keep.lib.sizes = FALSE] }
  design <- model.matrix(~group)
  v <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(v, design))
  res <- topTable(fit, coef = 2, number = Inf, sort.by = "P") %>%
    rownames_to_column("GeneID")
  write_csv(res, file.path(results_path, "limma", sprintf("limma_%s_%s.csv", patient_id, disease)))
}

# === Run all DE methods ===
for (i in seq_len(nrow(patients))) {
  patient_id <- patients$Patient_ID[i]
  disease <- patients$Disease[i]
  samples <- metadata %>% filter(Patient_ID == patient_id | Disease == "Healthy")
  
  run_deseq2(patient_id, disease, samples)
  run_edger(patient_id, disease, samples)
  run_limma(patient_id, disease, samples)
}

# === Merge DE results (significant genes only, shared across methods) ===
merge_and_filter <- function(file_path) {
  filename <- basename(file_path)
  patient_id <- str_match(filename, "^DESeq_s_([^_]+)")[,2]
  disease <- str_match(filename, "^DESeq_s_[^_]+_(.*)\\.csv$")[,2]
  
  edge_file <- file.path(results_path, "EdgeR", str_replace(filename, "^DESeq_", "edgeR_"))
  limma_file <- file.path(results_path, "limma", str_replace(filename, "^DESeq_", "limma_"))
  
  deseq2 <- read_csv(file_path, show_col_types = FALSE) %>%
    select(GeneID = 1, logFC = log2FoldChange, `P-value` = pvalue) %>%
    filter(`P-value` < 0.05) %>%
    mutate(Method = "DESeq2")
  
  edger <- read_csv(edge_file, show_col_types = FALSE) %>%
    select(GeneID, logFC, `P-value` = PValue) %>%
    filter(`P-value` < 0.05) %>%
    mutate(Method = "edgeR")
  
  limma <- read_csv(limma_file, show_col_types = FALSE) %>%
    select(GeneID, logFC, `P-value` = `P.Value`) %>%
    filter(`P-value` < 0.05) %>%
    mutate(Method = "limma_voom")
  
  common <- Reduce(intersect, list(deseq2$GeneID, edger$GeneID, limma$GeneID))
  
  bind_rows(deseq2, edger, limma) %>%
    filter(GeneID %in% common) %>%
    mutate(Patient = patient_id, Disease = disease) %>%
    select(GeneID, Patient, Disease, Method, logFC, `P-value`)
}

merged_results <- list.files(file.path(results_path, "DESeq2"), full.names = TRUE) %>%
  map_dfr(merge_and_filter)

write_csv(merged_results, file.path(results_path, "Merged", "DE_merged_filtered.csv"))
