library(limma)
library(edgeR)
library(dplyr)
library(readr)

# Paths
counts <- read.csv("D:/Research/GEO/GSE146853/filtered_counts.csv", check.names = FALSE, row.names = 1)
metadata <- read.csv("D:/Research/GEO/GSE146853/filtered_metadata.csv", check.names = FALSE, row.names = 1)
outdir <- "D:/Research/GEO/GSE146853/Results/limma"
dir.create(outdir, showWarnings = FALSE)
head(counts)
head(metadata)

# Loop through each patient (with replicates)
patients_to_test <- metadata %>%
  filter(Disease != "Healthy") %>%
  pull(Patient_ID) %>%
  unique()

for (patient_id in patients_to_test) {
  # Get disease name
  disease <- metadata %>%
    filter(Patient_ID == patient_id) %>%
    pull(Disease) %>%
    unique()
  
  # Select samples for this comparison
  selected_metadata <- metadata %>%
    filter(Patient_ID == patient_id | Disease == "Healthy")
  
  sample_ids <- selected_metadata$Sample_ID
  condition <- ifelse(selected_metadata$Disease == "Healthy", "Healthy", disease)
  group <- factor(condition, levels = c("Healthy", disease))
  
  # Subset count matrix
  counts_sub <- counts[, sample_ids]
  
  # Create DGEList
  dge <- DGEList(counts = counts_sub)
  dge <- calcNormFactors(dge)
  
  # Filter lowly expressed genes
  keep <- filterByExpr(dge, group = group)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # Design matrix
  design <- model.matrix(~ group)
  
  # Apply voom transformation
  v <- voom(dge, design, plot = FALSE)
  
  # Fit linear model
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  # Extract results: coef = 2 (disease vs Healthy)
  res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
  res <- tibble::rownames_to_column(res, var = "GeneID")
  
  # Save results
  out_file <- file.path(outdir, paste0("DEG_", patient_id, "_", disease, "_vs_Healthy_limma.csv"))
  write_csv(res, out_file)
  
  cat("limma+voom done:", patient_id, "-", disease, "\n")
}
