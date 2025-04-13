install.packages(c("tidyverse", "limma", "edgeR"))

library(tidyverse)
library(limma)
library(edgeR)

# Set file paths 
counts_file <- "D:/Research/GEO/GSE146853/raw_counts.csv"
metadata_file <- "D:/Research/GEO/GSE146853/metadata.csv"
output_counts_file <- "D:/Research/GEO/GSE146853/filtered_counts.csv"
output_metadata_file <- "D:/Research/GEO/GSE146853/filtered_metadata.csv"

# Metadata: clean and filter
metadata <- read.csv(metadata_file, row.names = "SampleName")
metadata$Sample_ID <- trimws(metadata$Sample_ID)
filtered_metadata <- metadata %>%
  mutate(Sample_ID = gsub("-", ".", gsub("_GeneCount", "", Sample_ID))) %>%
  filter(!is.na(Disease), Disease != "NA", Disease != "")

# Raw counts: clean and filter
raw_counts <- read.csv(counts_file, row.names = "GeneID")
colnames(raw_counts) <- trimws(colnames(raw_counts))
colnames(raw_counts) <- gsub("_GeneCount", "", colnames(raw_counts))
colnames(raw_counts) <- gsub("-", ".", colnames(raw_counts))

# Keep only samples that exist in both metadata and counts
common_samples <- intersect(colnames(raw_counts), filtered_metadata$Sample_ID)

# Check if any common samples exist
if (length(common_samples) == 0) {
  stop("No matching sample names found between raw_counts and metadata.")
}

raw_counts <- raw_counts[, common_samples]
filtered_metadata <- filtered_metadata %>% filter(Sample_ID %in% common_samples)
filtered_metadata <- filtered_metadata[match(colnames(raw_counts), filtered_metadata$Sample_ID), ]

# Create DGEList and filter lowly expressed genes
dge <- DGEList(counts = raw_counts, group = filtered_metadata$Disease)
keep <- filterByExpr(dge, group = dge$samples$group)
filtered_counts <- dge$counts[keep, ]

# Save filtered counts
write.csv(filtered_counts, file = output_counts_file)
write.csv(filtered_metadata, file = output_metadata_file)
