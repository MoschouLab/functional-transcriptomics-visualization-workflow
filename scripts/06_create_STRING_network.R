#!/usr/bin/env Rscript
################################################################################
# Script 06: STRING Protein-Protein Interaction Network
# Project: Nicotiana tabacum RNA-seq Day10 vs Day1 Analysis
# Author: Jose Moya Cuevas (Pepe) - MoschouLab
# Date: 2025-02-01
################################################################################

# Description:
# Creates protein-protein interaction network using STRING database
# Visualizes network with gene expression data overlay

################################################################################
# Load Required Libraries
################################################################################

suppressPackageStartupMessages({
  library(STRINGdb)
  library(igraph)
  library(ggraph)
  library(dplyr)
  library(tibble)
  library(readr)
  library(here)
})

cat("=======================================================================\n")
cat("SCRIPT 06: STRING PROTEIN-PROTEIN INTERACTION NETWORK\n")
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
n_top_genes <- 100  # Number of genes for network
string_species <- 3702  # Arabidopsis thaliana
confidence_score <- 150  # LOW confidence for more interactions (0-1000 scale)

# Color palette
color_up <- "#E74C3C"
color_down <- "#3498DB"

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

cat("\nStep 2: Selecting top genes for network...\n")

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

# Combine and clean Arabidopsis IDs
top_genes <- bind_rows(top_up, top_down) %>%
  filter(!is.na(at_accession) & at_accession != "") %>%
  mutate(
    # Remove isoform (.1, .2, etc.) - keep only locus (AT1G01010)
    at_locus = sub("\\.\\d+$", "", at_accession)
  ) %>%
  distinct(at_locus, .keep_all = TRUE)

cat(sprintf("  - Selected %d genes with Arabidopsis orthologs\n", nrow(top_genes)))
cat(sprintf("  - UP: %d, DOWN: %d\n", 
            sum(top_genes$regulation == 1), 
            sum(top_genes$regulation == 2)))

################################################################################
# Initialize STRING Database
################################################################################

cat("\nStep 3: Connecting to STRING database...\n")

string_db <- STRINGdb$new(
  version = "12.0",
  species = string_species,
  score_threshold = confidence_score,
  network_type = "full",
  input_directory = "~/stringdb_cache"
)

cat(sprintf("  - STRING database initialized (threshold: %d/1000)\n", confidence_score))

################################################################################
# Map Genes to STRING
################################################################################

cat("\nStep 4: Mapping genes to STRING identifiers...\n")

# Prepare data frame with cleaned locus IDs
genes_for_mapping <- data.frame(
  gene = top_genes$at_locus,
  stringsAsFactors = FALSE
)

# Map to STRING
mapped_genes <- string_db$map(
  genes_for_mapping,
  "gene",
  removeUnmappedRows = FALSE
)

# Remove unmapped genes
mapped_genes <- mapped_genes %>%
  filter(!is.na(STRING_id))

# Add back expression data
mapped_genes <- mapped_genes %>%
  left_join(
    top_genes %>% select(at_locus, gene_id, log2fc, regulation, regulation_label),
    by = c("gene" = "at_locus")
  )

cat(sprintf("  - Successfully mapped %d/%d genes to STRING (%.1f%%)\n", 
            nrow(mapped_genes), nrow(genes_for_mapping),
            100 * nrow(mapped_genes) / nrow(genes_for_mapping)))

# Check if we have enough genes
if (nrow(mapped_genes) < 10) {
  stop("Too few genes mapped to STRING (", nrow(mapped_genes), 
       "). Need at least 10 genes for network analysis.")
}

################################################################################
# Get Protein-Protein Interactions
################################################################################

cat("\nStep 5: Retrieving protein-protein interactions...\n")

# Get interactions
interactions <- string_db$get_interactions(mapped_genes$STRING_id)

if (nrow(interactions) == 0) {
  cat("  WARNING: No interactions found with threshold ", confidence_score, "/1000\n")
  cat("  Trying with lower threshold (100/1000)...\n")
  
  # Reinitialize with lower threshold
  string_db <- STRINGdb$new(
    version = "12.0",
    species = string_species,
    score_threshold = 100,
    network_type = "full",
    input_directory = "~/stringdb_cache"
  )
  
  interactions <- string_db$get_interactions(mapped_genes$STRING_id)
  
  if (nrow(interactions) == 0) {
    stop("No interactions found even with threshold 100/1000. Too few mapped genes.")
  }
  
  confidence_score <- 100
}

cat(sprintf("  - Retrieved %d interactions\n", nrow(interactions)))
cat(sprintf("  - Final confidence threshold: %d/1000\n", confidence_score))

################################################################################
# Build Network Graph
################################################################################

cat("\nStep 6: Building network graph...\n")

# Prepare edges with score
edges_df <- interactions %>%
  select(from, to, combined_score) %>%
  mutate(weight = combined_score / 1000)  # Normalize to 0-1

# Prepare vertices
vertices_df <- data.frame(
  name = mapped_genes$STRING_id,
  gene_id = mapped_genes$gene_id,
  log2fc = mapped_genes$log2fc,
  regulation = mapped_genes$regulation_label,
  stringsAsFactors = FALSE
)

# Create igraph object
g <- graph_from_data_frame(
  d = edges_df,
  directed = FALSE,
  vertices = vertices_df
)

# Remove self-loops and multiple edges
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

cat(sprintf("  - Nodes: %d\n", vcount(g)))
cat(sprintf("  - Edges: %d\n", ecount(g)))
cat(sprintf("  - Network density: %.3f\n", edge_density(g)))

# Calculate node metrics
V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g)

cat("  - Node metrics calculated\n")

################################################################################
# Create Network Visualization
################################################################################

cat("\nStep 7: Creating network visualization...\n")

# Prepare graph for ggraph
set.seed(123)  # For reproducibility

p <- ggraph(g, layout = "fr") +
  # Edges
  geom_edge_link(
    aes(width = weight),
    alpha = 0.3,
    color = "grey70"
  ) +
  scale_edge_width(range = c(0.3, 2), guide = "none") +
  # Nodes
  geom_node_point(
    aes(size = degree, color = log2fc),
    alpha = 0.8
  ) +
  scale_size_continuous(
    range = c(3, 12),
    name = "Degree\n(Connectivity)"
  ) +
  scale_color_gradient2(
    low = color_down,
    mid = "white",
    high = color_up,
    midpoint = 0,
    name = "log2FC"
  ) +
  # Labels for high-degree nodes
  geom_node_text(
    aes(filter = degree > quantile(degree, 0.85), 
        label = gene_id),
    size = 2.5,
    repel = TRUE,
    max.overlaps = 20
  ) +
  # Theme
  labs(
    title = "Protein-Protein Interaction Network",
    subtitle = sprintf("%d genes - STRING database (confidence >= %.2f)",
                      vcount(g), confidence_score/1000)
  ) +
  theme_graph(
    base_family = "sans",
    background = "white"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "right"
  )

# Save plot
output_file <- file.path(output_dir, "panel_F_STRING_network.pdf")

ggsave(
  output_file,
  plot = p,
  width = 14,
  height = 12,
  units = "in",
  dpi = 300
)

cat(sprintf("  Saved: %s\n", basename(output_file)))

################################################################################
# Network Statistics Summary
################################################################################

cat("\nStep 8: Network statistics summary...\n")

# Find hub genes (high degree)
hubs <- data.frame(
  gene_id = V(g)$gene_id,
  degree = V(g)$degree,
  log2fc = V(g)$log2fc,
  regulation = V(g)$regulation
) %>%
  arrange(desc(degree)) %>%
  head(10)

cat("\nTop 10 hub genes (highest connectivity):\n")
print(hubs)

# Community detection
communities <- cluster_louvain(g)
V(g)$community <- membership(communities)

cat(sprintf("\nNetwork communities detected: %d\n", length(communities)))

################################################################################
# Save Network Data
################################################################################

cat("\nStep 9: Saving network data...\n")

# Save node data
node_data <- data.frame(
  gene_id = V(g)$gene_id,
  STRING_id = V(g)$name,
  log2fc = V(g)$log2fc,
  regulation = V(g)$regulation,
  degree = V(g)$degree,
  betweenness = V(g)$betweenness,
  community = V(g)$community
) %>%
  arrange(desc(degree))

write_tsv(
  node_data,
  file.path(output_dir, "STRING_network_nodes.tsv")
)

# Save edge data
edge_data <- igraph::as_data_frame(g, "edges")

write_tsv(
  edge_data,
  file.path(output_dir, "STRING_network_edges.tsv")
)

cat("  - Saved node and edge data\n")

################################################################################
# Create Log File
################################################################################

log_file <- file.path(log_dir, paste0("06_STRING_network_", 
                                      format(Sys.time(), "%Y%m%d_%H%M%S"), 
                                      ".log"))
sink(log_file)
cat("STRING NETWORK ANALYSIS LOG\n")
cat("Execution time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Genes requested:", n_top_genes, "\n")
cat("Genes with orthologs:", nrow(top_genes), "\n")
cat("Mapped to STRING:", nrow(mapped_genes), "\n")
cat("Network nodes:", vcount(g), "\n")
cat("Network edges:", ecount(g), "\n")
cat("Network density:", edge_density(g), "\n")
cat("Confidence threshold:", confidence_score/1000, "\n")
cat("Communities detected:", length(communities), "\n\n")
cat("Top hub genes:\n")
print(hubs)
sink()

################################################################################
# Finish
################################################################################

cat("\n=======================================================================\n")
cat("SCRIPT 06 COMPLETED SUCCESSFULLY!\n")
cat(sprintf("Log saved: %s\n", basename(log_file)))
cat("=======================================================================\n\n")

cat("Files created:\n")
cat("  - panel_F_STRING_network.pdf (force-directed layout)\n")
cat("  - STRING_network_nodes.tsv (node data)\n")
cat("  - STRING_network_edges.tsv (edge data)\n\n")

cat("Next step: Assemble multipanel figure\n\n")
