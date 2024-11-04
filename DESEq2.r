#!/usr/bin/env Rscript

# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)

# Read count data and sample information
counts <- read.table("Count_matrix.txt", header=TRUE, row.names=1)
sample_info <- read.table("Sample_Info", header=TRUE, row.names=1)

# Create DESeqDataSet object with design formula
dds <- DESeqDataSetFromMatrix(countData=counts, colData=sample_info, design=~ tissue + time + tissue:time)

# Filter out low count genes
counts_filtered <- counts[rowSums(counts >= 5) > 0, ]

# Run DESeq2 analysis
dds <- DESeq(dds)

# Calculate size factors and normalized counts
sizeFactors(dds)
normalized_counts <- colSums(counts(dds, normalized=TRUE))

# Plot dispersion estimates
pdf("dispersion_estimates.pdf")  # Save plot as PDF
plotDispEsts(dds)
dev.off()  # Close the PDF device

# Get results and summarize
res <- results(dds)
summary(res)

# Order results by adjusted p-value
res <- res[order(res$padj),]

# Create a basic volcano plot
pdf("volcano_plot.pdf")  # Save plot as PDF
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(res, padj<.01), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
dev.off()  # Close the PDF device

# Variance-stabilizing transformation
vsdata <- vst(dds, blind=FALSE)

# PCA plot
pdf("PCA_plot.pdf")  # Save plot as PDF
plotPCA(vsdata, intgroup=c("tissue","time"))
dev.off()  # Close the PDF device

# Conduct contrasts for specific comparisons
res_tissue <- results(dds, contrast=c("tissue", "Liver", "Heart"))
res_time <- results(dds, contrast=c("time", "0_hrs", "12_hrs"))
res_interaction <- results(dds, name="tissueLiver.time12_hrs")

# Set thresholds for significance
padj_threshold <- 0.05
log2fc_threshold <- 2

# Filter significant genes for each contrast
res_tissue_sig <- as.data.frame(res_tissue) %>%
  filter(padj < padj_threshold & abs(log2FoldChange) > log2fc_threshold)

res_time_sig <- as.data.frame(res_time) %>%
  filter(padj < padj_threshold & abs(log2FoldChange) > log2fc_threshold)

res_interaction_sig <- as.data.frame(res_interaction) %>%
  filter(padj < padj_threshold & abs(log2FoldChange) > log2fc_threshold)

# Print the number of significant genes found
cat("Significant genes for tissue contrast: ", nrow(res_tissue_sig), "\n")
cat("Significant genes for time contrast: ", nrow(res_time_sig), "\n")
cat("Significant genes for interaction contrast: ", nrow(res_interaction_sig), "\n")

# Write results to CSV files
write.csv(res_tissue_sig, file = "significant_genes_tissue.csv", row.names = TRUE)
write.csv(res_time_sig, file = "significant_genes_time.csv", row.names = TRUE)
write.csv(res_interaction_sig, file = "significant_genes_interaction.csv", row.names = TRUE)

# Subset for specific comparisons
# Heart 0 hrs vs Heart 12 hrs
dds_heart <- dds[, dds$tissue == "Heart"]
res_Heart_0_vs_Heart_12 <- results(dds_heart, contrast=c("time", "12_hrs", "0_hrs"))

# Liver 0 hrs vs Liver 12 hrs
dds_liver <- dds[, dds$tissue == "Liver"]
res_Liver_0_vs_Liver_12 <- results(dds_liver, contrast=c("time", "12_hrs", "0_hrs"))

# Heart 0 hrs vs Liver 0 hrs
dds_0hrs <- dds[, dds$time == "0_hrs"]
res_Heart_0_vs_Liver_0 <- results(dds_0hrs, contrast=c("tissue", "Liver", "Heart"))

# Heart 12 hrs vs Liver 12 hrs
dds_12hrs <- dds[, dds$time == "12_hrs"]
res_Heart_12_vs_Liver_12 <- results(dds_12hrs, contrast=c("tissue", "Liver", "Heart"))

# Define a function to filter results based on thresholds
filter_results <- function(res) {
  res_df <- as.data.frame(res)
  res_filtered <- res_df[!is.na(res_df$padj) & res_df$padj < padj_threshold & 
                           (res_df$log2FoldChange > log2fc_threshold | res_df$log2FoldChange < -log2fc_threshold), ]
  return(res_filtered)
}

# Filter results for each comparison
res_Heart_0_vs_Heart_12_filtered <- filter_results(res_Heart_0_vs_Heart_12)
res_Liver_0_vs_Liver_12_filtered <- filter_results(res_Liver_0_vs_Liver_12)
res_Heart_0_vs_Liver_0_filtered <- filter_results(res_Heart_0_vs_Liver_0)
res_Heart_12_vs_Liver_12_filtered <- filter_results(res_Heart_12_vs_Liver_12)

# Save filtered results to CSV files
write.csv(res_Heart_0_vs_Heart_12_filtered, file="Filtered_DESeq2_Results_Heart_0_vs_Heart_12.csv", row.names=TRUE)
write.csv(res_Liver_0_vs_Liver_12_filtered, file="Filtered_DESeq2_Results_Liver_0_vs_Liver_12.csv", row.names=TRUE)
write.csv(res_Heart_0_vs_Liver_0_filtered, file="Filtered_DESeq2_Results_Heart_0_vs_Liver_0.csv", row.names=TRUE)
write.csv(res_Heart_12_vs_Liver_12_filtered, file="Filtered_DESeq2_Results_Heart_12_vs_Liver_12.csv", row.names=TRUE)

# Save filtered results to TSV files
write.table(res_Heart_0_vs_Heart_12_filtered, file="Filtered_DESeq2_Results_Heart_0_vs_Heart_12.tsv", sep="\t", row.names=TRUE, col.names=NA)
write.table(res_Liver_0_vs_Liver_12_filtered, file="Filtered_DESeq2_Results_Liver_0_vs_Liver_12.tsv", sep="\t", row.names=TRUE, col.names=NA)
write.table(res_Heart_0_vs_Liver_0_filtered, file="Filtered_DESeq2_Results_Heart_0_vs_Liver_0.tsv", sep="\t", row.names=TRUE, col.names=NA)
write.table(res_Heart_12_vs_Liver_12_filtered, file="Filtered_DESeq2_Results_Heart_12_vs_Liver_12.tsv", sep="\t", row.names=TRUE, col.names=NA)

# Extract top genes based on adjusted p-value
top_genes <- head(order(res$padj, na.last = NA), 50)

# Get normalized counts for top genes and replace NAs with 0
norm_counts <- counts(dds, normalized=TRUE)[top_genes, ]
norm_counts[is.na(norm_counts)] <- 0
norm_counts <- norm_counts[rowSums(norm_counts) > 0, ]

# Create heatmap of top differentially expressed genes
pdf("heatmap_top_genes.pdf")  # Save plot as PDF
pheatmap(norm_counts, 
         cluster_rows = TRUE,         # Cluster rows (genes)
         cluster_cols = TRUE,         # Cluster columns (samples)
         show_rownames = TRUE,        # Show gene names
         show_colnames = TRUE,        # Show sample names
         fontsize_row = 8,            # Font size for row names
         fontsize_col = 10,           # Font size for column names
         main = "Heatmap of Top 50 Differentially Expressed Genes", 
         scale = "row")
dev.off()  # Close the PDF device

