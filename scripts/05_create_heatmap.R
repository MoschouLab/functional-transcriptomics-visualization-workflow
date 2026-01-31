#!/usr/bin/env Rscript
################################################################################
# Script 05: Heatmap Creation with ComplexHeatmap
# Project: Nicotiana tabacum RNA-seq Day10 vs Day1 Analysis
# Author: Jose Moya Cuevas (Pepe) - MoschouLab
# Date: 2025-01-31
################################################################################

# Description:
# Creates publication-ready clustered heatmap of top DEGs

################################################################################
# Load Required Libraries
################################################################################

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tibble)
  library(readr)
  library(here)
})

cat("=======================================================================\n")
cat("SCRIPT 05: HEATMAP CREATION\n")
cat("=======================================================================\n\n")

################################################################################
# Define Paths and Parameters
################################################################################

input_dir <- here("data", "processed")
output_dir <- here("results", "figures")
log_dir <- here("logs")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

# Parameters
n_top_genes <- 100  # Total genes to show (50 UP + 50 DOWN)

# Color palette (matching previous plots)
color_up <- "#E74C3C"      # Red
color_down <- "#3498DB"    # Blue

################################################################################
# Load Data
################################################################################

cat("Step 1: Loading data...\n")

input_file <- file.path(input_dir, "degs_all_clean.tsv")
df <- read_tsv(input_file, show_col_types = FALSE)

cat(sprintf("  - Loaded %d DEGs\n", nrow(df)))

################################################################################
# Select Top Genes
################################################################################

cat("\nStep 2: Selecting top genes for heatmap...\n")

# Select top UP genes
top_up <- df %>%
  filter(regulation == 1) %>%
  arrange(desc(log2fc)) %>%
  head(n_top_genes / 2)

# Select top DOWN genes
top_down <- df %>%
  filter(regulation == 2) %>%
  arrange(log2fc) %>%
  head(n_top_genes / 2)

# Combine
top_genes <- bind_rows(top_up, top_down) %>%
  arrange(desc(log2fc))

cat(sprintf("  - Selected %d UP genes (log2FC: %.2f to %.2f)\n",
            nrow(top_up), min(top_up$log2fc), max(top_up$log2fc)))
cat(sprintf("  - Selected %d DOWN genes (log2FC: %.2f to %.2f)\n",
            nrow(top_down), min(top_down$log2fc), max(top_down$log2fc)))

################################################################################
# Prepare Expression Matrix
################################################################################

cat("\nStep 3: Preparing expression matrix...\n")

# Extract TPM columns
tpm_cols <- c("tpm_1d_rep1", "tpm_1d_rep2", "tpm_1d_rep3",
              "tpm_10d_rep1", "tpm_10d_rep2", "tpm_10d_rep3", "tpm_10d_rep4")

# Create matrix
expr_matrix <- top_genes %>%
  select(gene_id, all_of(tpm_cols)) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

# Log2 transform (add pseudocount to avoid log(0))
expr_matrix_log2 <- log2(expr_matrix + 1)

# Z-score normalization (by row)
expr_matrix_zscore <- t(scale(t(expr_matrix_log2)))

cat(sprintf("  - Expression matrix: %d genes x %d samples\n", 
            nrow(expr_matrix_zscore), ncol(expr_matrix_zscore)))
cat(sprintf("  - Z-score range: %.2f to %.2f\n",
            min(expr_matrix_zscore), max(expr_matrix_zscore)))

################################################################################
# Prepare Annotations
################################################################################

cat("\nStep 4: Preparing annotations...\n")

# Sample annotations (column)
sample_annotation <- data.frame(
  sample = colnames(expr_matrix_zscore),
  condition = c(rep("Day 1", 3), rep("Day 10", 4))
)

# Gene annotations (row)
gene_annotation <- top_genes %>%
  select(gene_id, regulation_label, log2fc) %>%
  column_to_rownames("gene_id")

cat("  - Sample annotation created\n")
cat("  - Gene annotation created\n")

################################################################################
# Create Heatmap
################################################################################

cat("\nStep 5: Creating compact heatmap with 100 genes...\n")

# Define color functions
col_fun <- colorRamp2(
  c(-3, 0, 3),
  c("#3498DB", "white", "#E74C3C")
)

# Column annotation (more compact)
col_ha <- HeatmapAnnotation(
  Condition = sample_annotation$condition,
  col = list(Condition = c("Day 1" = "#95A5A6", "Day 10" = "#2ECC71")),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10),
  annotation_height = unit(0.5, "cm"),
  annotation_legend_param = list(
    Condition = list(title = "Condition", title_gp = gpar(fontface = "bold", fontsize = 10))
  )
)

# Row annotation (more compact)
row_ha <- rowAnnotation(
  Regulation = gene_annotation$regulation_label,
  `log2FC` = anno_barplot(
    gene_annotation$log2fc,
    gp = gpar(fill = ifelse(gene_annotation$regulation_label == "UP", 
                            color_up, color_down)),
    width = unit(2, "cm"),
    axis_param = list(gp = gpar(fontsize = 8))
  ),
  col = list(Regulation = c("UP" = color_up, "DOWN" = color_down)),
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 10),
  annotation_width = unit(2.5, "cm"),
  annotation_legend_param = list(
    Regulation = list(title = "Regulation", title_gp = gpar(fontface = "bold", fontsize = 10))
  )
)

# Create heatmap - optimized for 100 genes in compact space
ht <- Heatmap(
  expr_matrix_zscore,
  name = "Z-score",
  col = col_fun,
  
  # Clustering
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "complete",
  
  # Show options - small font size for 100 genes
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),  # Very small but readable
  row_names_side = "left",
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 10),
  
  # Dendrogram - compact
  show_row_dend = TRUE,
  row_dend_width = unit(2, "cm"),
  
  # Annotations
  top_annotation = col_ha,
  right_annotation = row_ha,
  
  # Titles
  column_title = "Heatmap of Top 100 DEGs - Day 10 vs Day 1",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  
  # Legend
  heatmap_legend_param = list(
    title = "Z-score",
    title_gp = gpar(fontface = "bold", fontsize = 10),
    at = c(-3, -2, -1, 0, 1, 2, 3),
    labels = c("-3", "-2", "-1", "0", "1", "2", "3"),
    legend_height = unit(4, "cm")
  ),
  
  # Border
  border = TRUE,
  
  # Row split by regulation
  row_split = gene_annotation$regulation_label,
  row_title_gp = gpar(fontsize = 11, fontface = "bold"),
  
  # Minimal gap
  row_gap = unit(2, "mm")
)

# Save heatmap - SAME height as before (14 inches) but wider for readability
output_file <- file.path(output_dir, "panel_E_heatmap_top100.pdf")

cat("  Creating PDF with dimensions: 13 x 14 inches\n")
cat("  Note: 100 genes with font size 6pt - check if names are legible\n")

pdf(output_file, width = 13, height = 14)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

cat(sprintf("  Saved: %s\n", basename(output_file)))

################################################################################
# Summary Statistics
################################################################################

cat("\nStep 6: Creating summary statistics...\n")

summary_stats <- data.frame(
  metric = c("Total genes", "UP genes", "DOWN genes",
             "Mean log2FC (UP)", "Mean log2FC (DOWN)",
             "Max log2FC (UP)", "Min log2FC (DOWN)",
             "Z-score range", "Samples Day 1", "Samples Day 10"),
  value = c(nrow(top_genes), nrow(top_up), nrow(top_down),
            sprintf("%.2f", mean(top_up$log2fc)),
            sprintf("%.2f", mean(top_down$log2fc)),
            sprintf("%.2f", max(top_up$log2fc)),
            sprintf("%.2f", min(top_down$log2fc)),
            sprintf("%.2f to %.2f", min(expr_matrix_zscore), max(expr_matrix_zscore)),
            "3", "4")
)

cat("Heatmap summary:\n")
print(summary_stats)

################################################################################
# Create Log File
################################################################################

log_file <- file.path(log_dir, paste0("05_heatmap_", 
                                      format(Sys.time(), "%Y%m%d_%H%M%S"), 
                                      ".log"))
sink(log_file)
cat("HEATMAP CREATION LOG\n")
cat("Execution time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Number of genes:", n_top_genes, "\n")
cat("PDF dimensions: 13 x 14 inches\n")
cat("Font size: 6pt (compact layout)\n\n")
cat("Summary:\n")
print(summary_stats)
sink()

################################################################################
# Finish
################################################################################

cat("\n=======================================================================\n")
cat("SCRIPT 05 COMPLETED SUCCESSFULLY!\n")
cat(sprintf("Log saved: %s\n", basename(log_file)))
cat("=======================================================================\n\n")

cat("File created: panel_E_heatmap_top100.pdf (13 x 14 inches)\n")
cat("WARNING: Gene names use 6pt font - verify readability in PDF\n\n")

cat("Next step: Decide if Panel F (Network) should be created\n\n")
