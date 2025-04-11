install.packages(c("tidyverse", "limma", "edgeR"))

library(tidyverse)
library(limma)
library(edgeR)

# Set file paths
counts_file <- "D:/Research/GEO/GSE146853/raw_counts.csv"
metadata_file <- "D:/Research/GEO/GSE146853/metadata.csv"
output_file <- "D:/Research/GEO/GSE146853/filtered_counts.csv"

# Read input files
raw_counts <- read.csv(counts_file, row.names = 1)
metadata <- read.csv(metadata_file, sep = "\t")

# Sanitize column names (remove extra spaces, unify format)
colnames(raw_counts) <- trimws(colnames(raw_counts))
metadata$SampleName <- trimws(metadata$SampleName)

# Clean metadata (remove NA or invalid Disease samples)
metadata <- metadata %>%
  filter(!is.na(Disease), Disease != "NA", Disease != "")

# Keep only samples that exist in both metadata and counts
common_samples <- intersect(colnames(raw_counts), metadata$SampleName)

# Check if any common samples exist
if (length(common_samples) == 0) {
  stop("No matching sample names found between raw_counts and metadata.")
}

raw_counts <- raw_counts[, common_samples]
metadata <- metadata %>% filter(SampleName %in% common_samples)

# Ensure sample order matches
metadata <- metadata[match(colnames(raw_counts), metadata$SampleName), ]

# Create DGEList and filter lowly expressed genes
dge <- DGEList(counts = raw_counts, group = metadata$Disease)
keep <- filterByExpr(dge, group = dge$samples$group)
filtered_counts <- dge$counts[keep, ]

# Save filtered counts
write.csv(filtered_counts, file = output_file)
