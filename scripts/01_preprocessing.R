install.packages(c("tidyverse", "limma", "edgeR"))

library(tidyverse)
library(limma)
library(edgeR)
library(magrittr)

# Set file paths 
counts_file <- "D:/Research/GEO/GSE146853/raw_counts.csv"
metadata_file <- "D:/Research/GEO/GSE146853/metadata.csv"
output_counts_file <- "D:/Research/GEO/GSE146853/filtered_counts.csv"
output_metadata_file <- "D:/Research/GEO/GSE146853/filtered_metadata.csv"


# === Metadata ===
# Cleaning
metadata <- read.csv(metadata_file, row.names = "SampleName")
head(metadata)
metadata$Sample_ID <- trimws(metadata$Sample_ID)
metadata <- metadata %>%
  mutate(Sample_ID = gsub("-", ".", gsub("_GeneCount", "", Sample_ID))) %>%
  filter(!is.na(Disease), Disease != "NA", Disease != "")

# Filtering for Disease
metadata <- metadata %>%
  filter(!is.na(Disease), Disease != "NA", Disease != "")

# Add Patient ID column for replicate filtering
metadata$Patient_ID <- sub("\\.\\d+$", "", metadata$Sample_ID)

# Identify patients with exactly 2 replicates
replicated <- metadata %>%
  group_by(Patient_ID) %>%
  tally() %>%
  filter(n == 2) %>%
  pull(Patient_ID)

# Keep only replicated patients or healthy controls
filtered_metadata <- metadata %>%
  filter(Patient_ID %in% replicated)
head(filtered_metadata)


# === Raw counts ===
# Clean and filter
raw_counts <- read.csv(counts_file, row.names = "GeneID")
head(raw_counts)
colnames(raw_counts) <- colnames(raw_counts) %>%
  trimws() %>%
  gsub("_GeneCount", "", .) %>%
  gsub("-", ".", .)

# Keep only samples that exist in both metadata and counts
common_samples <- intersect(colnames(raw_counts), filtered_metadata$Sample_ID)

# Check if any common samples exist
if (length(common_samples) == 0) {
  stop("No matching sample names found between raw_counts and metadata.")
}

# Subset counts and metadata based on common samples
raw_counts <- raw_counts[, common_samples]
filtered_metadata <- filtered_metadata %>%
  filter(Sample_ID %in% common_samples) %>%
  arrange(match(Sample_ID, common_samples))  

# Check that everything matches
dim(filtered_metadata)
dim(raw_counts)


# === Create DGEList and filter lowly expressed genes ===
dge <- DGEList(counts = raw_counts, group = filtered_metadata$Disease)
keep <- filterByExpr(dge, group = dge$samples$group)
filtered_counts <- dge$counts[keep, ]

# Save filtered counts
write.csv(filtered_counts, file = output_counts_file)
write.csv(filtered_metadata, file = output_metadata_file)
