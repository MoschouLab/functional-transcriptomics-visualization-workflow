# Nicotiana tabacum RNA-seq Day10 vs Day1: Transcriptomics Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

## Overview

Comprehensive downstream transcriptomics analysis pipeline for *Nicotiana tabacum* RNA-seq data comparing Day 10 vs Day 1 developmental stages. This repository contains all scripts, data, and documentation for downstream bioinformatics analysis generating publication-ready figures.

## Project Structure
```
tobacco-rnaseq-day10-vs-day1/
├── data/
│   ├── raw/                    # Original DEG data
│   └── processed/              # Filtered DEGs by log2FC thresholds
├── results/
│   ├── figures/                # Publication-quality figures
│   ├── tables/                 # Enrichment analysis results
│   └── stats/                  # Summary statistics
├── scripts/                    # Analysis R scripts
├── logs/                       # Execution logs
├── environment.yml             # Conda environment specification
└── README.md                   # This file
```

## Dataset

- **Species**: *Nicotiana tabacum*
- **Comparison**: Day 10 vs Day 1
- **DEGs**: 1,389 differentially expressed genes
  - 675 UP-regulated (log2FC > 0)
  - 714 DOWN-regulated (log2FC < 0)
- **Significance**: FDR-adjusted p-value < 0.05

## Methodological Note: Use of Arabidopsis Orthologs

**Important**: *Nicotiana tabacum* lacks comprehensive GO/KEGG annotation databases in Bioconductor. Following standard practice in plant transcriptomics, this analysis uses **Arabidopsis thaliana orthologs** for functional enrichment analysis.

- **Ortholog mapping**: Pre-computed by RNA-seq pipeline (column "At_accesn" in input data)
- **Annotation source**: TAIR10 (The Arabidopsis Information Resource)
- **Justification**: 
  - *Arabidopsis* has the most complete plant functional annotation
  - High conservation of core biological processes across angiosperms
  - Bioconductor package `org.At.tair.db` provides comprehensive GO/KEGG mapping
- **Validation**: Ortholog assignments based on BLAST/OrthoFinder (upstream analysis)

This approach is widely accepted in plant omics studies for non-model species.

## Analyses Included

### Multi-panel Figure Components

- **Panel A**: GO enrichment circular plot (provided)
- **Panel B**: GO bubble plots comparison (BP/CC/MF) across log2FC thresholds (1.0, 1.5, 2.0)
- **Panel C**: KEGG pathway enrichment analysis
- **Panel D**: Volcano plot with functional annotation highlights
- **Panel E**: Heatmap of top differentially expressed genes
- **Panel F**: Protein-protein interaction network

## Installation

### 1. Clone the repository
```bash
git clone https://github.com/MoschouLab/tobacco-rnaseq-day10-vs-day1.git
cd tobacco-rnaseq-day10-vs-day1
```

### 2. Create conda environment
```bash
conda env create -f environment.yml
conda activate tobacco-rnaseq
```

### 3. Install Bioconductor packages
```R
# In R console
BiocManager::install(c(
  "clusterProfiler",
  "enrichplot", 
  "org.At.tair.db",
  "DOSE",
  "pathview",
  "ComplexHeatmap",
  "DESeq2"
))
```

## Usage

Run scripts sequentially:
```bash
# 1. Data preparation
Rscript scripts/01_data_preparation.R

# 2. GO enrichment analysis
Rscript scripts/02_GO_enrichment_analysis.R

# 3. KEGG pathway analysis
Rscript scripts/03_KEGG_pathway_analysis.R

# 4. Create volcano plot
Rscript scripts/04_create_volcano_plot.R

# 5. Create heatmap
Rscript scripts/05_create_heatmap.R

# 6. Network analysis
Rscript scripts/06_network_analysis.R

# 7. Assemble final figure
Rscript scripts/07_assemble_final_figure.R
```

## Results

All figures are saved in `results/figures/` as publication-ready PDFs.

## Citation

If you use this pipeline, please cite:
```
[To be added upon publication]
```

## FAIR Principles

This pipeline follows FAIR data principles:
- **Findable**: GitHub repository with DOI
- **Accessible**: Open-source MIT license
- **Interoperable**: Standard file formats (TSV, PDF, XLSX)
- **Reusable**: Documented code, conda environment, clear structure

## Contact

- **Supervisor**: Moschou PN
- **Postdoctoral Bioinformatician**: Moya-Cuevas J
- **Lab**: [MoschouLab](https://github.com/MoschouLab)
- **Institution**: IMBB-FORTH, University of Crete

## License

MIT License - see [LICENSE](LICENSE) file for details.
