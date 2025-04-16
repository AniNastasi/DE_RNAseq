if (!require(GGally)) install.packages("GGally")
if (!require(VennDiagram)) install.packages("VennDiagram")

# Load required libraries
library(tidyverse)
library(VennDiagram)
library(reshape2)
library(ggpubr)
library(GGally)
library(VennDiagram)

# Set path to your merged DE data
file_path <- "D:/Research/GEO/GSE146853/Results/Merged/DE_merged_filtered.csv"

# Load data
de_data <- read.csv(file_path)
head(de_data)

# Reshape data: wide format for method comparison
logfc_wide <- de_data %>%
  select(GeneID, Patient, Method, logFC) %>%
  pivot_wider(names_from = Method, values_from = logFC)

# Correlation plot between methods (logFC)
cor_plot <- ggpairs(logfc_wide[, c("DESeq2", "edgeR", "limma_voom")],
                    title = "Correlation of logFC between Methods")

# Show the plot
print(cor_plot)

# Determine significant genes (P < 0.05) per method
sig_genes <- de_data %>%
  filter(P.value < 0.05) %>%
  group_by(Method) %>%
  summarise(Significant_Genes = list(unique(GeneID)))

# Create Venn diagram
venn_list <- setNames(sig_genes$Significant_Genes, sig_genes$Method)

venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  fill = c("skyblue", "pink1", "lightgreen"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  main = "Significant Genes Overlap"
)

grid.newpage()
grid.draw(venn.plot)

# Summary stats: number of significant genes
sig_counts <- de_data %>%
  filter(P.value < 0.05) %>%
  group_by(Method) %>%
  summarise(Num_Significant = n_distinct(GeneID))

print(sig_counts)

# Optional: save outputs
ggsave("logFC_method_correlation_plot.png", plot = cor_plot, width = 8, height = 6)
# Venn can't be saved directly via ggsave; use Cairo or png
png("venn_significant_genes.png", width = 800, height = 800)
grid.draw(venn.plot)
dev.off()
