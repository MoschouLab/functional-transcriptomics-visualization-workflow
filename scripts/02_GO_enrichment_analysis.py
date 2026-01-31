#!/usr/bin/env python3
################################################################################
# Script 02: GO Enrichment Analysis using g:Profiler
# Project: Nicotiana tabacum RNA-seq Day10 vs Day1 Analysis
# Author: Jose Moya Cuevas (Pepe) - MoschouLab
# Date: 2025-01-31
################################################################################

"""
Performs GO enrichment analysis (BP, CC, MF) using g:Profiler API
Creates publication-ready bubble plots for UP and DOWN regulated genes
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from gprofiler import GProfiler
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

print("="*75)
print("SCRIPT 02: GO ENRICHMENT ANALYSIS (g:Profiler)")
print("="*75)
print()

################################################################################
# Configuration
################################################################################

# Paths
BASE_DIR = Path.cwd()
DATA_DIR = BASE_DIR / "data" / "processed"
RESULTS_FIGURES = BASE_DIR / "results" / "figures"
RESULTS_TABLES = BASE_DIR / "results" / "tables"
LOG_DIR = BASE_DIR / "logs"

RESULTS_FIGURES.mkdir(exist_ok=True, parents=True)
RESULTS_TABLES.mkdir(exist_ok=True, parents=True)
LOG_DIR.mkdir(exist_ok=True, parents=True)

# Use only one threshold (all genes)
THRESHOLD = 4.0

# GO categories
GO_CATEGORIES = {
    'BP': 'GO:BP',
    'CC': 'GO:CC',
    'MF': 'GO:MF'
}

GO_NAMES = {
    'BP': 'Biological Process',
    'CC': 'Cellular Component',
    'MF': 'Molecular Function'
}

# g:Profiler parameters
ORGANISM = 'athaliana'
SIGNIFICANCE_THRESHOLD = 0.05

# Color palette (matching original GO circular plot)
COLORS = {
    'UP': '#E74C3C',      # Red for UP-regulated
    'DOWN': '#3498DB',    # Blue for DOWN-regulated
}

# Plot parameters
TOP_N_TERMS = 20  # Number of top terms to show

################################################################################
# Load Data
################################################################################

print("Step 1: Loading data...")

file_path = DATA_DIR / f"degs_filtered_log2fc{THRESHOLD}.tsv"
df = pd.read_csv(file_path, sep='\t')
print(f"  - Loaded {len(df)} genes with |log2FC| >= {THRESHOLD}")

################################################################################
# Prepare Gene Lists
################################################################################

print("\nStep 2: Preparing gene lists for enrichment...")

# Filter genes with Arabidopsis orthologs
df_with_ortho = df[df['at_accession'].notna()].copy()

# UP-regulated genes
genes_up = df_with_ortho[df_with_ortho['regulation'] == 1]['at_accession'].tolist()

# DOWN-regulated genes
genes_down = df_with_ortho[df_with_ortho['regulation'] == 2]['at_accession'].tolist()

print(f"  UP-regulated: {len(genes_up)} genes")
print(f"  DOWN-regulated: {len(genes_down)} genes")
print()

################################################################################
# Perform g:Profiler Enrichment
################################################################################

print("Step 3: Performing g:Profiler enrichment analysis...")
print("This may take a few minutes...")
print()

gp = GProfiler(return_dataframe=True)

enrichment_results = {}

for ont_key, ont_source in GO_CATEGORIES.items():
    print(f"--- {GO_NAMES[ont_key]} ({ont_key}) ---")
    
    enrichment_results[ont_key] = {}
    
    # UP-regulated
    try:
        result_up = gp.profile(
            organism=ORGANISM,
            query=genes_up,
            sources=[ont_source],
            significance_threshold_method='fdr',
            user_threshold=SIGNIFICANCE_THRESHOLD,
            background=None,
        )
        
        if result_up is not None and len(result_up) > 0:
            enrichment_results[ont_key]['UP'] = result_up
            print(f"  UP: {len(result_up)} terms enriched")
        else:
            print("  UP: No enrichment")
            enrichment_results[ont_key]['UP'] = pd.DataFrame()
    except Exception as e:
        print(f"  UP: Error - {str(e)}")
        enrichment_results[ont_key]['UP'] = pd.DataFrame()
    
    # DOWN-regulated
    try:
        result_down = gp.profile(
            organism=ORGANISM,
            query=genes_down,
            sources=[ont_source],
            significance_threshold_method='fdr',
            user_threshold=SIGNIFICANCE_THRESHOLD,
            background=None,
        )
        
        if result_down is not None and len(result_down) > 0:
            enrichment_results[ont_key]['DOWN'] = result_down
            print(f"  DOWN: {len(result_down)} terms enriched")
        else:
            print("  DOWN: No enrichment")
            enrichment_results[ont_key]['DOWN'] = pd.DataFrame()
    except Exception as e:
        print(f"  DOWN: Error - {str(e)}")
        enrichment_results[ont_key]['DOWN'] = pd.DataFrame()
    
    print()

################################################################################
# Save Enrichment Results to Excel
################################################################################

print("Step 4: Saving enrichment results...")

for ont_key in GO_CATEGORIES.keys():
    output_file = RESULTS_TABLES / f"GO_enrichment_{ont_key}_results.xlsx"
    
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        for direction in ['UP', 'DOWN']:
            if direction in enrichment_results[ont_key] and len(enrichment_results[ont_key][direction]) > 0:
                df_result = enrichment_results[ont_key][direction]
                sheet_name = f"{direction}_regulated"
                df_result.to_excel(writer, sheet_name=sheet_name, index=False)
    
    print(f"  Saved: {output_file.name}")

print()

################################################################################
# Create Publication-Ready Bubble Plots (Separate for UP and DOWN)
################################################################################

print("Step 5: Creating publication-ready bubble plots...")

def create_single_bubble_plot(ont_key, direction):
    """
    Create individual bubble plot for UP or DOWN regulated genes
    """
    # Check if we have data
    if direction not in enrichment_results[ont_key] or len(enrichment_results[ont_key][direction]) == 0:
        print(f"    No enrichment for {direction}")
        return None
    
    # Get data
    df = enrichment_results[ont_key][direction].copy()
    
    # Calculate fold enrichment
    df['fold_enrichment'] = (df['intersection_size'] / df['query_size']) / \
                            (df['term_size'] / df['effective_domain_size'])
    
    df['neg_log10_fdr'] = -np.log10(df['p_value'])
    
    # Select top terms by p-value
    df_plot = df.nsmallest(TOP_N_TERMS, 'p_value').copy()
    
    # Truncate long term names
    df_plot['name_short'] = df_plot['name'].apply(
        lambda x: x[:70] + '...' if len(x) > 70 else x
    )
    
    # Sort by fold enrichment for better visualization
    df_plot = df_plot.sort_values('fold_enrichment')
    
    # Determine appropriate x-axis limits
    max_fe = df_plot['fold_enrichment'].max()
    min_fe = df_plot['fold_enrichment'].min()
    
    # Use log scale if fold enrichment range is very large
    use_log_scale = (max_fe / min_fe) > 50
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Create scatter plot
    scatter = ax.scatter(
        x=df_plot['fold_enrichment'],
        y=df_plot['name_short'],
        s=df_plot['intersection_size'] * 10,
        c=df_plot['neg_log10_fdr'],
        alpha=0.7,
        cmap='Reds' if direction == 'UP' else 'Blues',
        edgecolors='black',
        linewidth=0.8,
        vmin=2,
        vmax=df_plot['neg_log10_fdr'].max()
    )
    
    # Set x-scale
    if use_log_scale:
        ax.set_xscale('log')
        ax.set_xlabel('Fold Enrichment (log scale)', fontsize=12, fontweight='bold')
    else:
        ax.set_xlabel('Fold Enrichment', fontsize=12, fontweight='bold')
        ax.set_xlim(left=0, right=max_fe * 1.1)
    
    # Vertical reference line at fold enrichment = 1
    if not use_log_scale:
        ax.axvline(x=1, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    
    # Labels and title
    ax.set_ylabel('GO Term', fontsize=12, fontweight='bold')
    ax.set_title(f'{GO_NAMES[ont_key]} - {direction}-regulated Genes', 
                 fontsize=14, fontweight='bold', pad=15)
    
    # Grid
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5, axis='x')
    ax.set_axisbelow(True)
    
    # Tick parameters
    ax.tick_params(axis='both', labelsize=10)
    
    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label('-log10(FDR)', fontsize=11, fontweight='bold')
    cbar.ax.tick_params(labelsize=10)
    
    # Size legend
    from matplotlib.lines import Line2D
    legend_sizes = [10, 30, 50]
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', 
               markerfacecolor='gray', markersize=np.sqrt(s*10/np.pi), 
               alpha=0.7, markeredgecolor='black', markeredgewidth=0.8)
        for s in legend_sizes
    ]
    legend = ax.legend(
        legend_elements, 
        [str(s) for s in legend_sizes],
        title='Gene Count', 
        loc='lower right', 
        frameon=True, 
        fontsize=9,
        title_fontsize=10,
        handletextpad=1,
        borderpad=1,
        labelspacing=1.2
    )
    legend.get_title().set_fontweight('bold')
    
    plt.tight_layout()
    
    return fig

# Generate plots - separate files for UP and DOWN
for ont_key in GO_CATEGORIES.keys():
    print(f"  Creating plots for {GO_NAMES[ont_key]}...")
    
    for direction in ['UP', 'DOWN']:
        fig = create_single_bubble_plot(ont_key, direction)
        
        if fig is not None:
            output_file = RESULTS_FIGURES / f"panel_B_{ont_key}_{direction.lower()}_bubble.pdf"
            fig.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"    Saved: {output_file.name}")

print()

################################################################################
# Create Log File
################################################################################

from datetime import datetime

log_file = LOG_DIR / f"02_GO_enrichment_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

with open(log_file, 'w') as f:
    f.write("GO ENRICHMENT ANALYSIS LOG (g:Profiler)\n")
    f.write(f"Execution time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Organism: {ORGANISM} (Arabidopsis thaliana)\n")
    f.write(f"Threshold: |log2FC| >= {THRESHOLD}\n")
    f.write(f"GO categories: {', '.join(GO_CATEGORIES.keys())}\n")
    f.write(f"Significance threshold: {SIGNIFICANCE_THRESHOLD}\n")
    f.write(f"Top terms displayed: {TOP_N_TERMS}\n\n")
    
    f.write("Enrichment summary:\n")
    for ont_key in GO_CATEGORIES.keys():
        f.write(f"\n{GO_NAMES[ont_key]}:\n")
        for direction in ['UP', 'DOWN']:
            if direction in enrichment_results[ont_key] and len(enrichment_results[ont_key][direction]) > 0:
                n_terms = len(enrichment_results[ont_key][direction])
                f.write(f"  {direction}: {n_terms} terms\n")

print(f"Log saved: {log_file.name}")

################################################################################
# Finish
################################################################################

print()
print("="*75)
print("SCRIPT 02 COMPLETED SUCCESSFULLY!")
print("="*75)
print()
print("Files created:")
print(f"  - 3 Excel files with enrichment results (results/tables/)")
print(f"  - 6 bubble plot PDFs - UP and DOWN for each GO category (results/figures/)")
print()
print("Next step: Run script 03_KEGG_pathway_analysis.py")
print()
