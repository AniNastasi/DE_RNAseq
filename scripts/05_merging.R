library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)

# Define input/output paths
base_path <- "D:/Research/GEO/GSE146853/Results"
deseq2_path <- file.path(base_path, "DESeq2")
edger_path <- file.path(base_path, "EdgeR")
limma_path <- file.path(base_path, "limma")
output_file <- file.path(base_path, "Merged", "All_DEG_TidyFormat.csv")
dir.create(dirname(output_file), showWarnings = FALSE)

# List DESeq2 result files
deseq2_files <- list.files(deseq2_path, pattern = "^DEG_s_.*_vs_Healthy\\.csv$", full.names = TRUE)

# Function to read and tidy one patient's data
read_and_filter_tidy <- function(file_path) {
  filename <- basename(file_path)
  patient_id <- str_match(filename, "^DEG_s_([^_]+)")[,2]
  disease <- str_match(filename, "^DEG_s_[^_]+_([^_]+)_vs_Healthy")[,2]
  
  # Build corresponding file paths for edgeR and limma
  edger_file <- file.path(edger_path, paste0("DEG_s_", patient_id, "_", disease, "_vs_Healthy_edgeR.csv"))
  limma_file <- file.path(limma_path, paste0("DEG_s_", patient_id, "_", disease, "_vs_Healthy_limma.csv"))
  
  # === DESeq2 ===
  df_deseq2 <- suppressMessages(read_csv(file_path)) %>%
    select(GeneID = 1, logFC = log2FoldChange, `P-value` = pvalue) %>%
    filter(`P-value`< 0.05) %>%
    mutate(Method = "DESeq2")
  
  # === edgeR ===
  df_edger <- suppressMessages(read_csv(edger_file)) %>%
    select(GeneID, logFC = logFC, `P-value` = PValue) %>%
    filter(`P-value`< 0.05) %>%
    mutate(Method = "edgeR")
  
  # === limma+voom ===
  df_limma <- suppressMessages(read_csv(limma_file)) %>%
    select(GeneID, logFC = logFC, `P-value` = `P.Value`) %>%
    filter(`P-value`< 0.05) %>%
    mutate(Method = "limma_voom")
  
  # Find common genes across all 3 methods after filtering
  common_genes <- reduce(
    list(df_deseq2$GeneID, df_edger$GeneID, df_limma$GeneID),
    intersect
  )
  
  # Filter all to keep only shared significant genes
  df_all <- bind_rows(df_deseq2, df_edger, df_limma) %>%
    filter(GeneID %in% common_genes) %>%
    mutate(Patient = patient_id, Disease = disease) %>%
    select(GeneID, Patient, Disease, Method, logFC, `P-value`)
  
  return(df_all)
}

# Merge all into one long tidy dataframe
all_filtered_results <- map_dfr(deseq2_files, read_and_filter_tidy)

# Write to CSV
write_csv(all_filtered_results, output_file)

cat("✅ All DE results saved in tidy format to:\n", output_file, "\n")
