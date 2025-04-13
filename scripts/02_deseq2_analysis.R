
library(DESeq2)
library(dplyr)

counts <- read.csv("D:/Research/GEO/GSE146853/filtered_counts.csv", check.names = FALSE, row.names = 1)
metadata <- read.csv("D:/Research/GEO/GSE146853/filtered_metadata.csv", check.names = FALSE)
colnames(metadata) <- make.names(colnames(metadata), unique = TRUE)
metadata <- metadata[, c("Disease", "Sample_ID")]
metadata$Patient_ID <- sub("\\.\\d+$", "", metadata$Sample_ID)

# Identify patients with exactly 2 replicates
replicated <- metadata %>%
  group_by(Patient_ID) %>%
  tally() %>%
  filter(n == 2) %>%
  pull(Patient_ID)

# Filter metadata
metadata_filtered <- metadata %>%
  filter(Patient_ID %in% replicated | Disease == "Healthy")

# Output folder
outdir <- "D:/Research/GEO/GSE146853/Results/DESeq2"
dir.create(outdir, showWarnings = FALSE)

# Loop through each patient (with replicates)
patients_to_test <- metadata_filtered %>%
  filter(Disease != "Healthy") %>%
  pull(Patient_ID) %>%
  unique()

for (patient_id in patients_to_test) {
  disease <- metadata_filtered %>%
    filter(Patient_ID == patient_id) %>%
    pull(Disease) %>%
    unique()
  
  # Select samples: 2 for patient + all healthy
  selected_samples <- metadata_filtered %>%
    filter(Disease == "Healthy" | Patient_ID == patient_id)
  
  sample_ids <- selected_samples$Sample_ID
  cond <- ifelse(selected_samples$Disease == "Healthy", "Healthy", disease)
  condition <- factor(cond, levels = c("Healthy", disease))
  
  # Build colData
  coldata <- data.frame(row.names = sample_ids,
                        Condition = condition)
  
  # Subset count matrix
  counts_sub <- counts[, sample_ids]
  
  # Run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = round(counts_sub),
                                colData = coldata,
                                design = ~ Condition)
  dds <- DESeq(dds)
  res <- results(dds)
  
  # Save results
  out_file <- file.path(outdir, paste0("DEG_", patient_id, "_", disease, "_vs_Healthy.csv"))
  write.csv(as.data.frame(res), out_file)
  
  cat("DESeq2 done:", patient_id, "-", disease, "\n")
}
