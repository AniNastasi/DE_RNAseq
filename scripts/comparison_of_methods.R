required_packages <- c("tidyverse", "VennDiagram", "GGally", "gridExtra", "ggpubr", "reshape2")

new_packages <- required_packages[!(required_packages %in% installed.packages())]
if(length(new_packages)) install.packages(new_packages)
lapply(required_packages, library, character.only = TRUE)

# Set path and load data
file_path <- "D:/Research/GEO/GSE146853/Results/samples_with_replications/Merged/DE_merged_filtered.csv"
pictures_path <- "D:/Research/GEO/GSE146853/Results/samples_with_replications/pictures"
if (!dir.exists(pictures_path)) {
  dir.create(pictures_path, recursive = TRUE)
}
de_data <- read.csv(file_path)
head(de_data)

# Correlation of logFG between methods
logfc_wide <- de_data %>%
  select(GeneID, Patient, Method, logFC) %>%
  pivot_wider(names_from = Method, values_from = logFC)

cor_plot <- ggpairs(logfc_wide[, c("DESeq2", "edgeR", "limma_voom")],
                    title = "Correlation of logFC between Methods")
print(cor_plot)

# Spearman correlation by P-value ranks
avg_pvals <- de_data %>%
  group_by(Method, GeneID) %>%
  summarise(P_value = mean(P.value, na.rm = TRUE), .groups = "drop")
head(avg_pvals)

ranked <- avg_pvals %>%
  group_by(Method) %>%
  mutate(P_rank = rank(P_value, ties.method = "min")) %>%
  select(GeneID, Method, P_rank) %>%
  pivot_wider(names_from = Method, values_from = P_rank)

rank_cor <- cor(ranked[, -1], method = "spearman", use = "complete.obs")
print("Spearman Rank Correlation of Gene Significance:")
print(rank_cor)

# Jaccard similarity of top 100 significant genes per method
top_genes <- de_data %>%
  arrange(P.value) %>%
  group_by(Method) %>%
  slice_head(n = 100) %>%
  summarise(TopGenes = list(unique(GeneID)))

jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

jaccard_matrix <- matrix(NA, 3, 3)
methods <- top_genes$Method
for (i in 1:3) {
  for (j in 1:3) {
    jaccard_matrix[i, j] <- jaccard(top_genes$TopGenes[[i]], top_genes$TopGenes[[j]])
  }
}
rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- methods
print("Jaccard Similarity of Top 100 Genes:")
print(jaccard_matrix)

# Create Venn diagram of significantly DE genes
sig_genes <- de_data %>%
  filter(P.value < 0.05) %>%
  group_by(Method) %>%
  summarise(Significant_Genes = list(unique(GeneID)))

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

sig_counts <- de_data %>%
  filter(P.value < 0.05) %>%
  group_by(Method) %>%
  summarise(Num_Significant = n_distinct(GeneID))
print(sig_counts)

# Save pictures and matrices
ggsave(filename = file.path(pictures_path, "logFC_method_correlation_plot.png"),
       plot = cor_plot, width = 8, height = 6)

png(filename = file.path(pictures_path, "venn_significant_genes.png"),
    width = 800, height = 800)
grid.draw(venn.plot)
dev.off()

write.table(jaccard_matrix,
            file = file.path(pictures_path, "Jaccard_similarity_matrix.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

write.table(rank_cor,
            file = file.path(pictures_path, "Spearman_rank_correlation.txt"),
            sep = "\t", quote = FALSE, col.names = NA)