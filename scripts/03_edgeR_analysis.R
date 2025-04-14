library(edgeR)
library(dplyr)
library(readr)

# Paths
counts <- read.csv("D:/Research/GEO/GSE146853/filtered_counts.csv", check.names = FALSE, row.names = 1)
metadata <- read.csv("D:/Research/GEO/GSE146853/filtered_metadata.csv", check.names = FALSE, row.names = 1)
outdir <- "D:/Research/GEO/GSE146853/Results/EdgeR"
dir.create(outdir, showWarnings = FALSE)
head(counts)
head(metadata)

# Loop through each patient (with replicates)
patients_to_test <- metadata %>%
  filter(Disease != "Healthy") %>%
  pull(Patient_ID) %>%
  unique()
patients_to_test

for (patient_id in patients_to_test) {
  # Get disease label
  disease <- metadata %>%
    filter(Patient_ID == patient_id) %>%
    pull(Disease) %>%
    unique()
  
  # Select 2 patient replicates + all healthy samples
  selected_metadata <- metadata %>%
    filter(Patient_ID == patient_id | Disease == "Healthy")
  
  sample_ids <- selected_metadata$Sample_ID
  condition <- ifelse(selected_metadata$Disease == "Healthy", "Healthy", disease)
  group <- factor(condition, levels = c("Healthy", disease))
  
  # Subset count data
  counts_sub <- counts[, sample_ids]
  
  # edgeR: create DGEList
  dge <- DGEList(counts = counts_sub, group = group)
  dge <- calcNormFactors(dge)
  
  # Filter lowly expressed genes (again, just in case)
  keep <- filterByExpr(dge, group = dge$samples$group)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # Design matrix and model
  design <- model.matrix(~ group)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = 2)  # coef = 2: disease vs Healthy
  
  # Extract and save top DEGs
  results <- topTags(qlf, n = Inf)$table
  results <- tibble::rownames_to_column(results, var = "GeneID")
  
  # Output file
  out_file <- file.path(outdir, paste0("DEG_", patient_id, "_", disease, "_vs_Healthy_edgeR.csv"))
  write_csv(results, out_file)
  
  cat("edgeR done:", patient_id, "-", disease, "\n")
}
