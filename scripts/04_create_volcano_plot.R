#!/usr/bin/env Rscript
################################################################################
# Script 04: Volcano Plot Creation
# Project: Nicotiana tabacum RNA-seq Day10 vs Day1 Analysis
# Author: Jose Moya Cuevas (Pepe) - MoschouLab
# Date: 2025-01-31
################################################################################

# Description:
# Creates publication-ready volcano plot with annotated top genes

################################################################################
# Load Required Libraries
################################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(readr)
  library(here)
})

cat("=======================================================================\n")
cat("SCRIPT 04: VOLCANO PLOT CREATION\n")
cat("=======================================================================\n\n")

################################################################################
# Define Paths and Parameters
################################################################################

input_dir <- here("data", "processed")
output_dir <- here("results", "figures")
log_dir <- here("logs")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

# Thresholds
log2fc_threshold <- 4.0
qval_threshold <- 0.05

# Color palette (matching previous plots)
color_up <- "#E74C3C"      # Red
color_down <- "#3498DB"    # Blue
color_ns <- "#95A5A6"      # Grey

# Number of genes to label
n_top_genes <- 15  # Top genes by |log2FC| to annotate

################################################################################
# Load Data
################################################################################

cat("Step 1: Loading data...\n")

input_file <- file.path(input_dir, "degs_all_clean.tsv")
df <- read_tsv(input_file, show_col_types = FALSE)

cat(sprintf("  - Loaded %d DEGs\n", nrow(df)))
cat(sprintf("  - log2FC range: %.2f to %.2f\n", min(df$log2fc), max(df$log2fc)))
cat(sprintf("  - qvalue range: %.2e to %.2e\n\n", min(df$qvalue), max(df$qvalue)))

################################################################################
# Prepare Data for Plotting
################################################################################

cat("Step 2: Preparing data for volcano plot...\n")

# Calculate -log10(qvalue)
df <- df %>%
  mutate(
    neg_log10_qval = -log10(qvalue),
    # Classification
    regulation_status = case_when(
      regulation == 1 ~ "UP",
      regulation == 2 ~ "DOWN",
      TRUE ~ "NS"
    ),
    # For labeling: select top genes by |log2FC|
    abs_log2fc = abs(log2fc)
  )

# Select top genes to label
top_up <- df %>%
  filter(regulation == 1) %>%
  arrange(desc(log2fc)) %>%
  head(n_top_genes) %>%
  pull(gene_id)

top_down <- df %>%
  filter(regulation == 2) %>%
  arrange(log2fc) %>%
  head(n_top_genes) %>%
  pull(gene_id)

genes_to_label <- c(top_up, top_down)

df <- df %>%
  mutate(
    to_label = gene_id %in% genes_to_label
  )

cat(sprintf("  - UP-regulated: %d genes\n", sum(df$regulation == 1)))
cat(sprintf("  - DOWN-regulated: %d genes\n", sum(df$regulation == 2)))
cat(sprintf("  - Genes to annotate: %d\n\n", length(genes_to_label)))

################################################################################
# Create Volcano Plot
################################################################################

cat("Step 3: Creating volcano plot...\n")

# Calculate appropriate x-axis limits
x_min <- floor(min(df$log2fc)) - 0.5
x_max <- ceiling(max(df$log2fc)) + 0.5

cat(sprintf("  - X-axis range: %.1f to %.1f\n", x_min, x_max))

p <- ggplot(df, aes(x = log2fc, y = neg_log10_qval)) +
  # Points
  geom_point(
    aes(color = regulation_status),
    alpha = 0.6,
    size = 1.5
  ) +
  # Color scale
  scale_color_manual(
    values = c("UP" = color_up, "DOWN" = color_down, "NS" = color_ns),
    name = "Regulation"
  ) +
  # Set appropriate x-axis limits
  scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = seq(ceiling(x_min), floor(x_max), by = 2)
  ) +
  # Threshold lines
  geom_vline(
    xintercept = c(-log2fc_threshold, log2fc_threshold),
    linetype = "dashed",
    color = "black",
    linewidth = 0.5,
    alpha = 0.7
  ) +
  geom_hline(
    yintercept = -log10(qval_threshold),
    linetype = "dashed",
    color = "black",
    linewidth = 0.5,
    alpha = 0.7
  ) +
  # Gene labels
  geom_text_repel(
    data = df %>% filter(to_label),
    aes(label = gene_id),
    size = 2.5,
    max.overlaps = 30,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.size = 0.3,
    segment.color = "grey50",
    min.segment.length = 0,
    force = 2,
    force_pull = 0.5
  ) +
  # Labels and theme
  labs(
    title = "Volcano Plot - Day 10 vs Day 1",
    subtitle = sprintf("Threshold: |log2FC| >= %.1f, q-value < %.2f", 
                      log2fc_threshold, qval_threshold),
    x = expression(log[2]~"Fold Change"),
    y = expression(-log[10]~"q-value")
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5)
  )

# Save plot
output_file <- file.path(output_dir, "panel_D_volcano_plot.pdf")
ggsave(
  output_file,
  plot = p,
  width = 10,
  height = 10,
  units = "in",
  dpi = 300
)

cat(sprintf("  Saved: %s\n\n", basename(output_file)))

################################################################################
# Create Summary Statistics
################################################################################

cat("Step 4: Creating summary statistics...\n")

summary_stats <- df %>%
  group_by(regulation_status) %>%
  summarise(
    n_genes = n(),
    mean_log2fc = mean(log2fc),
    median_log2fc = median(log2fc),
    min_log2fc = min(log2fc),
    max_log2fc = max(log2fc),
    mean_qval = mean(qvalue),
    median_qval = median(qvalue),
    .groups = "drop"
  )

cat("Summary statistics by regulation status:\n")
print(summary_stats)

################################################################################
# Create Log File
################################################################################

log_file <- file.path(log_dir, paste0("04_volcano_plot_", 
                                      format(Sys.time(), "%Y%m%d_%H%M%S"), 
                                      ".log"))
sink(log_file)
cat("VOLCANO PLOT CREATION LOG\n")
cat("Execution time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Thresholds: |log2FC| >=", log2fc_threshold, ", q-value <", qval_threshold, "\n")
cat("X-axis range:", x_min, "to", x_max, "\n")
cat("\nGene counts:\n")
print(summary_stats)
sink()

################################################################################
# Finish
################################################################################

cat("\n=======================================================================\n")
cat("SCRIPT 04 COMPLETED SUCCESSFULLY!\n")
cat(sprintf("Log saved: %s\n", basename(log_file)))
cat("=======================================================================\n\n")

cat("File created: panel_D_volcano_plot.pdf\n\n")

cat("Next step: Run script 05_create_heatmap.R\n\n")
