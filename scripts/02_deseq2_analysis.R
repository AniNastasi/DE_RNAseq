library(DESeq2)
library(tidyverse)

counts <- read.csv("D:/Research/GEO/GSE146853/filtered_counts.csv", row.names = 1)
metadata <- read.csv("D:/Research/GEO/GSE146853/metadata.csv", sep = "\t")
head(counts)
head(metadata)

patients <- metadata %>% filter(Disease != "Healthy" & Disease != "NA")
controls <- metadata %>% filter(Disease == "Healthy")

for (i in 1:nrow(patients)) {
  patient_sample <- patients$SampleName[i]
  patient_id <- gsub("[^A-Za-z0-9]", "_", patient_sample)
  subset_samples <- c(patient_sample, controls$SampleName)
  
  sub_counts <- counts[, subset_samples]
  sub_meta <- metadata %>% filter(SampleName %in% subset_samples) %>%
    mutate(group = ifelse(SampleName == patient_sample, "Patient", "Control"))
  
  rownames(sub_meta) <- sub_meta$SampleName
  sub_meta <- sub_meta[colnames(sub_counts), ]  # Fix ordering
  
  dds <- DESeqDataSetFromMatrix(countData = sub_counts, colData = sub_meta, design = ~ group)
  dds <- DESeq(dds)
  res <- results(dds)
  
  write.csv(as.data.frame(res), file = paste0("D:/Research/GEO/GSE146853/Results/deseq2_", patient_id, ".csv"))
}

cat("DESeq2 patient-wise comparisons complete.\n")