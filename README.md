# Functional Transcriptomics Visualization Workflow

<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) -->
<!-- TODO: Update with actual Zenodo DOI after first release -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python Version](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![R Version](https://img.shields.io/badge/R-4.3.0+-blue.svg)](https://www.r-project.org/)

A comprehensive workflow for visualization and functional analysis of differentially expressed genes (DEGs) from RNA-seq experiments. This repository provides reproducible scripts for generating publication-ready figures including GO/KEGG enrichment analysis, volcano plots, expression heatmaps, and protein-protein interaction networks.

---

## 📖 Overview

This workflow analyzes differentially expressed genes from *Nicotiana tabacum* (tobacco) RNA-seq data comparing Day 10 vs Day 1 samples. The pipeline generates six publication-ready figure panels:

- **Panel A**: Gene Ontology circular plot (original)
- **Panel B**: GO enrichment bubble plots (BP/CC/MF × UP/DOWN)
- **Panel C**: KEGG pathway diverging barplot
- **Panel D**: Volcano plot with top DEG annotations
- **Panel E**: Hierarchical clustering heatmap (Z-score normalized)
- **Panel F**: STRING protein-protein interaction network

---

## ✨ Features

- 🧬 **Functional Enrichment**: GO and KEGG pathway analysis using g:Profiler
- 📊 **Publication-Ready Visualizations**: High-quality PDF outputs at 300 DPI
- 🔗 **Network Analysis**: Protein-protein interactions from STRING database
- 🎨 **Colorblind-Friendly**: Accessible color schemes throughout
- 📈 **Flexible Thresholds**: Easy adjustment of fold-change and significance cutoffs
- 🔄 **Reproducible**: Conda/pip environments with pinned versions

---

## 📁 Repository Structure
```
.
├── scripts/                    # Analysis scripts (R and Python)
│   ├── 01_data_preparation.R
│   ├── 02_GO_enrichment_analysis.py
│   ├── 03_KEGG_pathway_analysis.py
│   ├── 04_create_volcano_plot.R
│   ├── 05_create_heatmap.R
│   └── 06_create_STRING_network.R
├── data/
│   ├── raw/                   # Input data (Excel with DEGs)
│   └── processed/             # Cleaned TSV files
├── results/
│   ├── figures/               # Output plots (PDF, gitignored)
│   ├── tables/                # Enrichment results (XLSX/TSV, gitignored)
│   └── stats/                 # Summary statistics
├── logs/                      # Execution logs (gitignored)
├── environment.yml            # Conda environment
├── requirements.txt           # Python dependencies
└── README.md
```

---

## 🚀 Installation

### Prerequisites

- **R** (≥ 4.3.0)
- **Python** (≥ 3.10)
- **Conda** (recommended) or pip

### Option 1: Conda Environment (Recommended)
```bash
# Clone the repository
git clone https://github.com/MoschouLab/functional-transcriptomics-visualization-workflow.git
cd functional-transcriptomics-visualization-workflow

# Create conda environment
conda env create -f environment.yml
conda activate tobacco-rnaseq

# Install R packages (if not already installed)
Rscript -e "install.packages(c('ggplot2', 'ggrepel', 'ComplexHeatmap', 'circlize', 'dplyr', 'readr', 'here', 'tibble'), repos='https://cloud.r-project.org')"
Rscript -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('STRINGdb')"
```

### Option 2: pip + Manual R Installation
```bash
# Python dependencies
pip install -r requirements.txt

# R packages (run in R console)
install.packages(c('ggplot2', 'ggrepel', 'ComplexHeatmap', 'circlize', 'dplyr', 'readr', 'here', 'tibble'))
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("STRINGdb")
```

---

## 💻 Usage

### Input Data Format

Place your DEG analysis results in `data/raw/DEGs_Day10_vs_Day1.xlsx` with the following columns:

- `gene_id`: Gene identifier
- `log2fc`: Log2 fold-change
- `pvalue`: P-value
- `qvalue`: Adjusted p-value (FDR)
- `tpm_*`: TPM values for each replicate
- `tair_description`: Gene description (optional)
- `at_accession`: Arabidopsis ortholog (for enrichment analysis)

### Running the Workflow

Execute scripts sequentially:
```bash
# 1. Data preparation (filter by thresholds)
Rscript scripts/01_data_preparation.R

# 2. GO enrichment analysis
python scripts/02_GO_enrichment_analysis.py

# 3. KEGG pathway analysis
python scripts/03_KEGG_pathway_analysis.py

# 4. Volcano plot
Rscript scripts/04_create_volcano_plot.R

# 5. Expression heatmap
Rscript scripts/05_create_heatmap.R

# 6. PPI network (requires STRING database download, ~5-15 min first time)
Rscript scripts/06_create_STRING_network.R
```

### Outputs

All figures are saved as PDF files in `results/figures/`:
- `panel_B_*.pdf` (6 GO bubble plots)
- `panel_C_KEGG_diverging_barplot.pdf`
- `panel_D_volcano_plot.pdf`
- `panel_E_heatmap_top50.pdf` and `panel_E_heatmap_top100.pdf`
- `panel_F_STRING_network.pdf`

Enrichment results are saved as Excel/TSV files in `results/tables/`.

---

## ⚙️ Key Parameters

### Script 01: Data Preparation
- **Thresholds**: `|log2FC| >= 4.0, 5.0, 6.0` (default: 4.0)
- **q-value**: `< 0.05`

### Script 02-03: Enrichment Analysis
- **Organism**: Arabidopsis thaliana (3702)
- **FDR threshold**: `< 0.05`
- **Background**: All genes (no custom background)

### Script 05: Heatmap
- **Top genes**: 50 or 100 (25/50 UP + 25/50 DOWN)
- **Normalization**: Z-score by gene
- **Clustering**: Euclidean distance, complete linkage

### Script 06: PPI Network
- **Confidence score**: `≥ 0.15` (low-medium confidence)
- **Hub threshold**: Top 15% by degree
- **Layout**: Fruchterman-Reingold force-directed

---

## 📦 Software Requirements

### Python Packages
- pandas (≥ 2.0.3)
- matplotlib (≥ 3.7.1)
- seaborn (≥ 0.12.2)
- gprofiler-official
- openpyxl

### R Packages
- ggplot2 (≥ 3.4.2)
- ggrepel (≥ 0.9.3)
- ComplexHeatmap (≥ 2.10.0)
- igraph (≥ 1.5.0)
- ggraph (≥ 2.1.0)
- STRINGdb (≥ 2.14.0)
- dplyr, readr, here, tibble, circlize

---

## 🛠️ Troubleshooting

### STRING Database Download Timeout

If Script 06 times out downloading the STRING database:
```bash
# Pre-download manually
mkdir -p ~/stringdb_cache
cd ~/stringdb_cache
wget https://stringdb-downloads.org/download/protein.links.v12.0/3702.protein.links.v12.0.txt.gz

# Update script to use cache (line 62)
# Change: input_directory = ""
# To:     input_directory = "~/stringdb_cache"
```

### Memory Issues with Heatmap

For large gene sets (>200), use the top 50 version instead of top 100:
```bash
# Generates panel_E_heatmap_top50.pdf only
# Comment out top 100 section in script 05
```

---

## 🤝 Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

### Reporting Issues

Please use [GitHub Issues](https://github.com/MoschouLab/functional-transcriptomics-visualization-workflow/issues) to report bugs or request features.

---

## 📄 Citation

If you use this workflow in your research, please cite:
```bibtex
@software{functional_transcriptomics_workflow,
  author = {Moya-Cuevas, J and Moschou, PN.},
  title = {Functional Transcriptomics Visualization Workflow},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/MoschouLab/functional-transcriptomics-visualization-workflow},
  doi = {10.5281/zenodo.XXXXXXX},
  note = {Version 1.0.0}
}
```

> **Note**: The Zenodo DOI badge and citation will be updated after the first GitHub release is synchronized with Zenodo. To obtain a DOI:
> 1. Create a GitHub release (v1.0.0)
> 2. Zenodo will automatically create an archive and assign a DOI
> 3. Update the badge at the top of this README
> 4. Update the `doi` field in the citation

**Related publication** (manuscript in preparation):

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

### Third-Party Data

- **g:Profiler**: [GPL-3.0](https://biit.cs.ut.ee/gprofiler/)
- **STRING**: [CC BY 4.0](https://string-db.org/cgi/access?footer_active_subpage=licensing)
- **GO Consortium**: [CC BY 4.0](http://geneontology.org/docs/go-citation-policy/)
- **KEGG**: [Academic use](https://www.kegg.jp/kegg/legal.html)

---

## 👥 Authors

**Lead Bioinformatics Scientist**:
- **Jose Moya-Cuevas** - Postdoctoral Researcher
  - IMBB-FORTH, University of Crete
  - ORCID: [0000-0001-9537-8556](https://orcid.org/0000-0001-9537-8556)
  - Email: jose_moya@imbb.forth.gr

**Principal Investigator**:
- **Panagiotis N. Moschou** - Associate Professor
  - IMBB-FORTH, University of Crete
  - ERC Consolidator Grant PLANTEX
  - ORCID: [0000-0001-7212-0595](https://orcid.org/0000-0001-7212-0595)
  - Email: panagiotis.moschou@imbb.forth.gr
  - Lab: [Moschou Lab](https://pmoschoulab.org)

---

## 🙏 Acknowledgments

- **ERC Consolidator Grant PLANTEX** for funding
- **g:Profiler**, **STRING**, **GO Consortium**, and **KEGG** for databases
- **Bioconductor** and **CRAN** communities for excellent R packages

---

## 📚 Related Resources

- [Moschou Lab Website](https://pmoschoulab.org/)
- [g:Profiler Documentation](https://biit.cs.ut.ee/gprofiler/page/docs)
- [STRING Database](https://string-db.org/)
- [ComplexHeatmap Book](https://jokergoo.github.io/ComplexHeatmap-reference/book/)

---

**Version**: 1.0.0  
**Last Updated**: 2026-02-01  
**Status**: Production-ready

---

**Keywords**: RNA-seq, transcriptomics, differential expression, functional enrichment, GO analysis, KEGG pathways, protein-protein interactions, data visualization, bioinformatics workflow, *Nicotiana tabacum*
