#!/usr/bin/env Rscript
################################################################################
# Script 01: Data Preparation and Filtering
# Project: Nicotiana tabacum RNA-seq Day10 vs Day1 Analysis
# Author: Jose Moya Cuevas (Pepe) - MoschouLab
# Date: 2025-01-31
################################################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(here)
})

cat("=======================================================================\n")
cat("SCRIPT 01: DATA PREPARATION\n")
cat("=======================================================================\n\n")

################################################################################
# Define Paths
################################################################################

input_file <- here("data", "raw", "DEGs_Day10_vs_Day1.xlsx")
output_dir_processed <- here("data", "processed")
output_dir_stats <- here("results", "stats")
log_dir <- here("logs")

dir.create(output_dir_processed, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir_stats, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

################################################################################
# Load Data
################################################################################

cat("Step 1: Loading raw data...\n")

raw_data <- read_excel(
  input_file,
  sheet = "Stat_day10_X_Exp_day1",
  skip = 1
)

cat(sprintf("  - Loaded %d genes with %d columns\n", nrow(raw_data), ncol(raw_data)))

################################################################################
# Clean and Rename Columns
################################################################################

cat("\nStep 2: Cleaning data...\n")

degs_clean <- raw_data %>%
  rename(
    gene_id = gene_ID,
    ncbi_refseq = `NCBI Ref Seq`,
    eff_length = eff_length,
    tpm_1d_rep1 = KM1921,
    tpm_1d_rep2 = KM1922,
    tpm_1d_rep3 = KM1923,
    tpm_10d_rep1 = KM1924,
    tpm_10d_rep2 = KM1918,
    tpm_10d_rep3 = KM1920,
    tpm_10d_rep4 = KM1928,
    mean_tpm_1d = `1d`,
    mean_tpm_10d = `10d`,
    pvalue = pval,
    qvalue = qval,
    log2fc = b,
    regulation = `10DO vs 1DO`,
    tair_description = `TAIR10_description,  Edwards annotation, or Solanum lycopersium note`,
    at_accession = At_accesn,
    at_evalue = At_evalue
  ) %>%
  select(-target_id, -Note, -starts_with("..."))

# Convert to appropriate types
degs_clean <- degs_clean %>%
  mutate(
    across(starts_with("tpm_"), as.numeric),
    across(c(mean_tpm_1d, mean_tpm_10d), as.numeric),
    pvalue = as.numeric(pvalue),
    qvalue = as.numeric(qvalue),
    log2fc = as.numeric(log2fc),
    regulation = as.integer(regulation),
    eff_length = as.numeric(eff_length)
  )

# Check for NA in critical columns
n_na_log2fc <- sum(is.na(degs_clean$log2fc))
if (n_na_log2fc > 0) {
  cat(sprintf("  WARNING: Found %d gene(s) with NA in log2fc column\n", n_na_log2fc))
  na_genes <- degs_clean %>% 
    filter(is.na(log2fc)) %>% 
    pull(gene_id)
  cat(sprintf("  Genes with NA log2fc: %s\n", paste(na_genes, collapse=", ")))
  cat("  Removing these genes from analysis...\n")
  
  degs_clean <- degs_clean %>% 
    filter(!is.na(log2fc))
}

# Add regulation labels
degs_clean <- degs_clean %>%
  mutate(
    regulation_label = case_when(
      regulation == 1 ~ "UP",
      regulation == 2 ~ "DOWN",
      TRUE ~ NA_character_
    )
  )

cat(sprintf("  - Clean dataset: %d genes\n", nrow(degs_clean)))

################################################################################
# Generate Summary Statistics
################################################################################

cat("\nStep 3: Generating summary statistics...\n")

summary_stats <- list(
  total_degs = nrow(degs_clean),
  up_regulated = sum(degs_clean$regulation == 1, na.rm = TRUE),
  down_regulated = sum(degs_clean$regulation == 2, na.rm = TRUE),
  log2fc_range = range(degs_clean$log2fc, na.rm = TRUE),
  log2fc_mean = mean(degs_clean$log2fc, na.rm = TRUE),
  log2fc_median = median(degs_clean$log2fc, na.rm = TRUE),
  abs_log2fc_min = min(abs(degs_clean$log2fc), na.rm = TRUE),
  pvalue_range = range(degs_clean$pvalue, na.rm = TRUE),
  qvalue_range = range(degs_clean$qvalue, na.rm = TRUE),
  with_arabidopsis_ortholog = sum(!is.na(degs_clean$at_accession))
)

cat("\n--- DATASET SUMMARY ---\n")
cat(sprintf("Total DEGs: %d\n", summary_stats$total_degs))
cat(sprintf("  UP-regulated: %d (%.1f%%)\n", 
            summary_stats$up_regulated,
            100 * summary_stats$up_regulated / summary_stats$total_degs))
cat(sprintf("  DOWN-regulated: %d (%.1f%%)\n", 
            summary_stats$down_regulated,
            100 * summary_stats$down_regulated / summary_stats$total_degs))
cat(sprintf("log2FC range: %.3f to %.3f\n", 
            summary_stats$log2fc_range[1], 
            summary_stats$log2fc_range[2]))
cat(sprintf("|log2FC| minimum: %.3f (all genes are highly DE)\n", 
            summary_stats$abs_log2fc_min))
cat(sprintf("log2FC mean: %.3f, median: %.3f\n", 
            summary_stats$log2fc_mean,
            summary_stats$log2fc_median))
cat(sprintf("Genes with Arabidopsis orthologs: %d (%.1f%%)\n",
            summary_stats$with_arabidopsis_ortholog,
            100 * summary_stats$with_arabidopsis_ortholog / summary_stats$total_degs))

################################################################################
# Filter by log2FC Thresholds
################################################################################

cat("\nStep 4: Creating filtered datasets by log2FC thresholds...\n")
cat("NOTE: Dataset is pre-filtered with high |log2FC| values.\n")
cat("      Using thresholds: 4.0, 5.0, 6.0 (appropriate for this data)\n\n")

# Define thresholds adjusted for this dataset
thresholds <- c(4.0, 5.0, 6.0)
filtered_datasets <- list()

for (threshold in thresholds) {
  
  filtered <- degs_clean %>%
    filter(abs(log2fc) >= threshold) %>%
    arrange(desc(abs(log2fc)))
  
  filtered_datasets[[as.character(threshold)]] <- filtered
  
  n_up <- sum(filtered$regulation == 1, na.rm = TRUE)
  n_down <- sum(filtered$regulation == 2, na.rm = TRUE)
  
  cat(sprintf("  |log2FC| >= %.1f: %d genes (%d UP, %d DOWN)\n",
              threshold, nrow(filtered), n_up, n_down))
  
  output_file <- file.path(
    output_dir_processed,
    sprintf("degs_filtered_log2fc%.1f.tsv", threshold)
  )
  
  write_tsv(filtered, output_file)
  cat(sprintf("    Saved: %s\n", basename(output_file)))
}

################################################################################
# Save Full Dataset
################################################################################

cat("\nStep 5: Saving complete cleaned dataset...\n")

full_output <- file.path(output_dir_processed, "degs_all_clean.tsv")
write_tsv(degs_clean, full_output)
cat(sprintf("  Saved: %s (%d genes)\n", basename(full_output), nrow(degs_clean)))

################################################################################
# Save Summary Statistics
################################################################################

cat("\nStep 6: Saving summary statistics...\n")

threshold_summary <- data.frame(
  threshold = c("All", paste0(">=", thresholds)),
  total_genes = c(
    nrow(degs_clean),
    sapply(filtered_datasets, nrow)
  ),
  up_regulated = c(
    summary_stats$up_regulated,
    sapply(filtered_datasets, function(df) sum(df$regulation == 1, na.rm = TRUE))
  ),
  down_regulated = c(
    summary_stats$down_regulated,
    sapply(filtered_datasets, function(df) sum(df$regulation == 2, na.rm = TRUE))
  )
) %>%
  mutate(
    pct_up = round(100 * up_regulated / total_genes, 1),
    pct_down = round(100 * down_regulated / total_genes, 1)
  )

stats_file <- file.path(output_dir_stats, "summary_statistics.tsv")
write_tsv(threshold_summary, stats_file)

stats_txt <- file.path(output_dir_stats, "summary_statistics.txt")
sink(stats_txt)
cat("================================================================================\n")
cat("NICOTIANA TABACUM RNA-SEQ ANALYSIS - SUMMARY STATISTICS\n")
cat("Comparison: Day 10 vs Day 1\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("================================================================================\n\n")

cat("OVERALL DATASET:\n")
cat(sprintf("  Total DEGs: %d\n", summary_stats$total_degs))
cat(sprintf("  UP-regulated (log2FC > 0): %d (%.1f%%)\n",
            summary_stats$up_regulated,
            100 * summary_stats$up_regulated / summary_stats$total_degs))
cat(sprintf("  DOWN-regulated (log2FC < 0): %d (%.1f%%)\n",
            summary_stats$down_regulated,
            100 * summary_stats$down_regulated / summary_stats$total_degs))
cat(sprintf("\nlog2FC Statistics:\n"))
cat(sprintf("  Range: %.3f to %.3f\n", summary_stats$log2fc_range[1], summary_stats$log2fc_range[2]))
cat(sprintf("  Mean: %.3f\n", summary_stats$log2fc_mean))
cat(sprintf("  Median: %.3f\n", summary_stats$log2fc_median))
cat(sprintf("  Minimum |log2FC|: %.3f\n", summary_stats$abs_log2fc_min))
cat("\nNOTE: Dataset is pre-filtered - all genes have |log2FC| >= 4.0\n")
cat("      This represents very strong differential expression\n")

cat("\n--------------------------------------------------------------------------------\n")
cat("FILTERED DATASETS BY |log2FC| THRESHOLD:\n")
cat("(Thresholds adjusted to 4.0, 5.0, 6.0 based on data distribution)\n")
cat("--------------------------------------------------------------------------------\n\n")

print(threshold_summary)

cat("\n================================================================================\n")
cat("ANNOTATION COVERAGE:\n")
cat("================================================================================\n")
cat(sprintf("Genes with Arabidopsis orthologs (At_accession): %d / %d (%.1f%%)\n",
            summary_stats$with_arabidopsis_ortholog,
            summary_stats$total_degs,
            100 * summary_stats$with_arabidopsis_ortholog / summary_stats$total_degs))

cat("\n================================================================================\n")
sink()

cat(sprintf("  Saved: %s\n", basename(stats_file)))
cat(sprintf("  Saved: %s\n", basename(stats_txt)))

################################################################################
# Create Log File
################################################################################

log_file <- file.path(log_dir, paste0("01_data_preparation_", 
                                      format(Sys.time(), "%Y%m%d_%H%M%S"), 
                                      ".log"))
sink(log_file)
cat("DATA PREPARATION LOG\n")
cat("Execution time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total genes processed:", nrow(degs_clean), "\n")
cat(sprintf("Genes removed (NA in log2fc): %d\n", n_na_log2fc))
cat("Thresholds used: 4.0, 5.0, 6.0 (adjusted for pre-filtered data)\n")
cat("\nFiles created:\n")
cat("  - data/processed/degs_all_clean.tsv\n")
for (threshold in thresholds) {
  cat(sprintf("  - data/processed/degs_filtered_log2fc%.1f.tsv\n", threshold))
}
cat("  - results/stats/summary_statistics.tsv\n")
cat("  - results/stats/summary_statistics.txt\n")
sink()

cat("\n=======================================================================\n")
cat("SCRIPT 01 COMPLETED SUCCESSFULLY!\n")
cat(sprintf("Log saved: %s\n", basename(log_file)))
cat("=======================================================================\n\n")
cat("Next step: Run script 02_GO_enrichment_analysis.R\n\n")
