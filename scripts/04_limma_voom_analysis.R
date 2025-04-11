library(limma)
library(edgeR)
library(tidyverse)

counts <- read.csv("D:/Research/STW5/GEO/Results/filtered_counts.csv", row.names = 1)
metadata <- read.csv("D:/Research/STW5/GEO/metadata.csv", sep = "\t")

patients <- metadata %>% filter(Disease != "Healthy" & Disease != "NA")
controls <- metadata %>% filter(Disease == "Healthy")

for (i in 1:nrow(patients)) {
  patient_sample <- patients$SampleName[i]
  patient_id <- gsub("[^A-Za-z0-9]", "_", patient_sample)
  subset_samples <- c(patient_sample, controls$SampleName)
  
  sub_counts <- counts[, subset_samples]
  sub_meta <- metadata %>% filter(SampleName %in% subset_samples) %>%
    mutate(group = ifelse(SampleName == patient_sample, "Patient", "Control"))
  
  dge <- DGEList(counts = sub_counts)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~ group, data = sub_meta)
  v <- voom(dge, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = 2, number = Inf)
  
  write.csv(res, file = paste0("D:/Research/STW5/GEO/Results/limma_voom_", patient_id, ".csv"))
}

cat("limma+voom patient-wise comparisons complete.\n")