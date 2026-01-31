# Analysis Log - Nicotiana tabacum RNA-seq Day10 vs Day1

## Project Information

- **Project Name**: Tobacco RNA-seq Transcriptomics Analysis
- **Analysis Date**: January 31, 2025
- **Postdoctoral Bioinformatician**: Jose Moya Cuevas
- **Lab**: MoschouLab, University of Crete
- **Organism**: *Nicotiana tabacum*
- **Comparison**: Day 10 vs Day 1

---

## Dataset Summary

### Input Data
- **Source File**: `Copy_of_RNAseq_tobacco_cultures_10vs1_DEGs.xlsx`
- **Total DEGs**: 1,389 genes
- **UP-regulated**: 675 genes (log2FC > 0)
- **DOWN-regulated**: 714 genes (log2FC < 0)
- **Significance cutoff**: FDR-adjusted p-value < 0.05

### Sample Information
- **Condition 1 (Day1)**: 3 replicates (KM1921, KM1922, KM1923)
- **Condition 2 (Day10)**: 3 replicates (KM1924, KM1918, KM1920)
- **Expression metric**: TPM (Transcripts Per Million)

---

## Analysis Pipeline

### Phase 1: Project Setup (Jan 31, 2025)

**Completed:**
- ✅ Created project directory structure following FAIR principles
- ✅ Initialized git repository
- ✅ Created conda environment `tobacco-rnaseq`
- ✅ Established FAIR compliance files:
  - environment.yml
  - requirements.txt
  - .gitignore
  - LICENSE (MIT)
  - README.md
  - GitHub Actions CI workflow
- ✅ Data imported and validated

**Environment specifications:**
- Python 3.10
- R 4.3
- Key packages: tidyverse, clusterProfiler, ggplot2, pheatmap

---

### Phase 2: Data Preparation (In Progress)

**Goal**: Process raw DEG data and create filtered datasets for different log2FC thresholds

**Tasks:**
- [ ] Script 01: Data preparation
  - [ ] Load and clean Excel data
  - [ ] Create filtered datasets (log2FC: 1.0, 1.5, 2.0)
  - [ ] Generate summary statistics
  - [ ] Export processed data to TSV format

---

### Phase 3: Enrichment Analysis (Planned)

**Tasks:**
- [ ] Script 02: GO enrichment analysis
  - [ ] BP, CC, MF categories
  - [ ] Multiple log2FC thresholds comparison
  - [ ] Generate bubble plots
  
- [ ] Script 03: KEGG pathway analysis
  - [ ] Pathway enrichment
  - [ ] UP vs DOWN comparison
  - [ ] Generate visualization

---

### Phase 4: Figure Generation (Planned)

**Multi-panel figure components:**
- [ ] Panel A: GO circular plot (provided by collaborator)
- [ ] Panel B: GO bubble plots (3 thresholds × 3 categories)
- [ ] Panel C: KEGG pathway enrichment
- [ ] Panel D: Annotated volcano plot
- [ ] Panel E: Heatmap of top DEGs
- [ ] Panel F: PPI network analysis
- [ ] Final: Assembled multi-panel figure

---

## Technical Notes

### Color Palette (from original GO circular plot)
- **UP-regulated genes**: Red/Pink gradient (#E74C3C, #FF6B6B)
- **DOWN-regulated genes**: Blue/Purple gradient (#3498DB, #9B59B6)
- **Z-score gradient**: Light pink (increasing) → Purple (decreasing)

### Analysis Considerations
- Single comparison (Day10 vs Day1) - no multi-condition Venn diagrams
- Focus on threshold-dependent enrichment changes
- Emphasis on biological process interpretation
- Network analysis for functional relationships

---

## Issues and Solutions

### Issue 1: Excel file structure
- **Problem**: Multi-header Excel format
- **Solution**: Manual inspection and proper header parsing in R

---

## Next Steps

1. Install Bioconductor packages in R
2. Develop Script 01 for data preparation
3. Validate processed data structure
4. Proceed with enrichment analysis

---

## References

- Original GO plot and input dataset provided by collaborator
- Following FAIR principles

---

*Last updated: January 31, 2025*
