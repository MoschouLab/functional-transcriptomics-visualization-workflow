#!/usr/bin/env python3
################################################################################
# Script 03: KEGG Pathway Enrichment Analysis using g:Profiler
# Project: Nicotiana tabacum RNA-seq Day10 vs Day1 Analysis
# Author: Jose Moya Cuevas (Pepe) - MoschouLab
# Date: 2025-01-31
################################################################################

"""
Performs KEGG pathway enrichment analysis using g:Profiler API
Creates publication-ready diverging barplot for UP and DOWN regulated genes
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
print("SCRIPT 03: KEGG PATHWAY ENRICHMENT ANALYSIS (g:Profiler)")
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

# g:Profiler parameters
ORGANISM = 'athaliana'
SIGNIFICANCE_THRESHOLD = 0.05

# Color palette (matching original GO circular plot)
COLOR_UP = '#E74C3C'      # Red for UP-regulated
COLOR_DOWN = '#3498DB'    # Blue for DOWN-regulated

# Plot parameters
MAX_PATHWAYS_PER_DIRECTION = 15  # Maximum pathways to show per direction

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
# Perform g:Profiler KEGG Enrichment
################################################################################

print("Step 3: Performing KEGG pathway enrichment analysis...")
print("This may take a few minutes...")
print()

gp = GProfiler(return_dataframe=True)

enrichment_results = {}

# UP-regulated
try:
    result_up = gp.profile(
        organism=ORGANISM,
        query=genes_up,
        sources=['KEGG'],
        significance_threshold_method='fdr',
        user_threshold=SIGNIFICANCE_THRESHOLD,
        background=None,
    )
    
    if result_up is not None and len(result_up) > 0:
        enrichment_results['UP'] = result_up
        print(f"  UP-regulated genes: {len(result_up)} pathways enriched")
    else:
        print("  UP-regulated genes: No enrichment")
        enrichment_results['UP'] = pd.DataFrame()
except Exception as e:
    print(f"  UP-regulated genes: Error - {str(e)}")
    enrichment_results['UP'] = pd.DataFrame()

# DOWN-regulated
try:
    result_down = gp.profile(
        organism=ORGANISM,
        query=genes_down,
        sources=['KEGG'],
        significance_threshold_method='fdr',
        user_threshold=SIGNIFICANCE_THRESHOLD,
        background=None,
    )
    
    if result_down is not None and len(result_down) > 0:
        enrichment_results['DOWN'] = result_down
        print(f"  DOWN-regulated genes: {len(result_down)} pathways enriched")
    else:
        print("  DOWN-regulated genes: No enrichment")
        enrichment_results['DOWN'] = pd.DataFrame()
except Exception as e:
    print(f"  DOWN-regulated genes: Error - {str(e)}")
    enrichment_results['DOWN'] = pd.DataFrame()

print()

################################################################################
# Save Enrichment Results to Excel
################################################################################

print("Step 4: Saving enrichment results...")

output_file = RESULTS_TABLES / "KEGG_pathway_enrichment_results.xlsx"

with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
    for direction in ['UP', 'DOWN']:
        if direction in enrichment_results and len(enrichment_results[direction]) > 0:
            df_result = enrichment_results[direction]
            sheet_name = f"{direction}_regulated"
            df_result.to_excel(writer, sheet_name=sheet_name, index=False)

print(f"  Saved: {output_file.name}")
print()

################################################################################
# Create Diverging Barplot
################################################################################

print("Step 5: Creating publication-ready diverging barplot...")

def create_diverging_barplot():
    """
    Create a diverging barplot showing UP (right) and DOWN (left) pathways
    """
    # Check if we have any data
    has_up = 'UP' in enrichment_results and len(enrichment_results['UP']) > 0
    has_down = 'DOWN' in enrichment_results and len(enrichment_results['DOWN']) > 0
    
    if not has_up and not has_down:
        print("  No enrichment data to plot")
        return None
    
    # Prepare data
    plot_data = []
    
    # Process UP-regulated pathways
    if has_up:
        df_up = enrichment_results['UP'].copy()
        
        # Calculate fold enrichment
        df_up['fold_enrichment'] = (df_up['intersection_size'] / df_up['query_size']) / \
                                    (df_up['term_size'] / df_up['effective_domain_size'])
        
        df_up['neg_log10_fdr'] = -np.log10(df_up['p_value'])
        
        # Select top pathways
        df_up = df_up.nsmallest(min(MAX_PATHWAYS_PER_DIRECTION, len(df_up)), 'p_value')
        
        # Clean names
        df_up['pathway'] = df_up['name'].str.replace('Arabidopsis thaliana:', '').str.strip()
        df_up['pathway_short'] = df_up['pathway'].apply(
            lambda x: x[:60] + '...' if len(x) > 60 else x
        )
        
        df_up['direction'] = 'UP'
        plot_data.append(df_up[['pathway_short', 'fold_enrichment', 'intersection_size', 
                                 'neg_log10_fdr', 'direction']])
    
    # Process DOWN-regulated pathways
    if has_down:
        df_down = enrichment_results['DOWN'].copy()
        
        # Calculate fold enrichment (negative for left side)
        df_down['fold_enrichment'] = -((df_down['intersection_size'] / df_down['query_size']) / \
                                       (df_down['term_size'] / df_down['effective_domain_size']))
        
        df_down['neg_log10_fdr'] = -np.log10(df_down['p_value'])
        
        # Select top pathways
        df_down = df_down.nsmallest(min(MAX_PATHWAYS_PER_DIRECTION, len(df_down)), 'p_value')
        
        # Clean names
        df_down['pathway'] = df_down['name'].str.replace('Arabidopsis thaliana:', '').str.strip()
        df_down['pathway_short'] = df_down['pathway'].apply(
            lambda x: x[:60] + '...' if len(x) > 60 else x
        )
        
        df_down['direction'] = 'DOWN'
        plot_data.append(df_down[['pathway_short', 'fold_enrichment', 'intersection_size', 
                                   'neg_log10_fdr', 'direction']])
    
    # Combine data
    df_plot = pd.concat(plot_data, ignore_index=True)
    
    # Sort by fold enrichment (so DOWN is at bottom, UP at top)
    df_plot = df_plot.sort_values('fold_enrichment', ascending=True)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, max(6, len(df_plot) * 0.4)))
    
    # Create barplot
    y_pos = np.arange(len(df_plot))
    
    # Color by direction and intensity by FDR
    colors = []
    for _, row in df_plot.iterrows():
        if row['direction'] == 'UP':
            # Red gradient based on significance
            alpha = min(1.0, row['neg_log10_fdr'] / df_plot['neg_log10_fdr'].max())
            colors.append((0.91, 0.30, 0.24, 0.5 + 0.5*alpha))  # Red with alpha
        else:
            # Blue gradient based on significance
            alpha = min(1.0, row['neg_log10_fdr'] / df_plot['neg_log10_fdr'].max())
            colors.append((0.20, 0.60, 0.86, 0.5 + 0.5*alpha))  # Blue with alpha
    
    bars = ax.barh(y_pos, df_plot['fold_enrichment'], color=colors, 
                   edgecolor='black', linewidth=0.8)
    
    # Add gene count labels on bars
    for i, (idx, row) in enumerate(df_plot.iterrows()):
        x_pos = row['fold_enrichment']
        gene_count = int(row['intersection_size'])
        
        # Position label inside or outside bar depending on bar length
        if abs(x_pos) > 2:
            # Inside bar
            x_text = x_pos - (0.3 if x_pos > 0 else -0.3)
            color = 'white'
        else:
            # Outside bar
            x_text = x_pos + (0.3 if x_pos > 0 else -0.3)
            color = 'black'
        
        ax.text(x_text, i, f'{gene_count}', 
               va='center', ha='center', 
               fontsize=9, fontweight='bold', color=color)
    
    # Customize axes
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df_plot['pathway_short'], fontsize=10)
    ax.set_xlabel('Fold Enrichment', fontsize=12, fontweight='bold')
    ax.set_ylabel('KEGG Pathway', fontsize=12, fontweight='bold')
    
    # Add vertical line at x=0
    ax.axvline(x=0, color='black', linewidth=1.5, linestyle='-')
    
    # Add labels for UP/DOWN
    y_max = len(df_plot)
    if has_up:
        ax.text(ax.get_xlim()[1] * 0.85, y_max * 0.95, 'UP-regulated', 
               fontsize=11, fontweight='bold', color=COLOR_UP,
               bbox=dict(boxstyle='round,pad=0.5', facecolor='white', 
                        edgecolor=COLOR_UP, linewidth=2))
    
    if has_down:
        ax.text(ax.get_xlim()[0] * 0.85, y_max * 0.05, 'DOWN-regulated', 
               fontsize=11, fontweight='bold', color=COLOR_DOWN,
               bbox=dict(boxstyle='round,pad=0.5', facecolor='white', 
                        edgecolor=COLOR_DOWN, linewidth=2))
    
    # Title
    ax.set_title('KEGG Pathway Enrichment - UP vs DOWN Regulated Genes', 
                 fontsize=14, fontweight='bold', pad=20)
    
    # Grid
    ax.grid(True, axis='x', alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)
    
    # Spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    
    return fig

# Create diverging barplot
fig = create_diverging_barplot()

if fig is not None:
    output_file = RESULTS_FIGURES / "panel_C_KEGG_diverging_barplot.pdf"
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {output_file.name}")
else:
    print("  No plot generated (no enrichment data)")

print()

################################################################################
# Create Alternative: Individual Barplots
################################################################################

print("Step 6: Creating individual barplots for UP and DOWN...")

def create_individual_barplot(direction):
    """
    Create individual barplot for UP or DOWN
    """
    if direction not in enrichment_results or len(enrichment_results[direction]) == 0:
        print(f"    No enrichment for {direction}")
        return None
    
    df = enrichment_results[direction].copy()
    
    # Calculate metrics
    df['fold_enrichment'] = (df['intersection_size'] / df['query_size']) / \
                            (df['term_size'] / df['effective_domain_size'])
    df['neg_log10_fdr'] = -np.log10(df['p_value'])
    
    # Select top
    df = df.nsmallest(min(MAX_PATHWAYS_PER_DIRECTION, len(df)), 'p_value')
    
    # Clean names
    df['pathway'] = df['name'].str.replace('Arabidopsis thaliana:', '').str.strip()
    df['pathway_short'] = df['pathway'].apply(
        lambda x: x[:65] + '...' if len(x) > 65 else x
    )
    
    # Sort by fold enrichment
    df = df.sort_values('fold_enrichment', ascending=True)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, max(4, len(df) * 0.45)))
    
    y_pos = np.arange(len(df))
    
    # Color gradient by significance
    color_base = COLOR_UP if direction == 'UP' else COLOR_DOWN
    colors = []
    for fdr in df['neg_log10_fdr']:
        alpha = min(1.0, fdr / df['neg_log10_fdr'].max())
        if direction == 'UP':
            colors.append((0.91, 0.30, 0.24, 0.5 + 0.5*alpha))
        else:
            colors.append((0.20, 0.60, 0.86, 0.5 + 0.5*alpha))
    
    bars = ax.barh(y_pos, df['fold_enrichment'], color=colors,
                   edgecolor='black', linewidth=0.8)
    
    # Add gene count labels
    for i, (idx, row) in enumerate(df.iterrows()):
        gene_count = int(row['intersection_size'])
        x_text = row['fold_enrichment'] + 0.3
        ax.text(x_text, i, f'{gene_count}', 
               va='center', ha='left', fontsize=9, fontweight='bold')
    
    # Customize
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df['pathway_short'], fontsize=10)
    ax.set_xlabel('Fold Enrichment', fontsize=11, fontweight='bold')
    ax.set_ylabel('KEGG Pathway', fontsize=11, fontweight='bold')
    ax.set_title(f'KEGG Pathway Enrichment - {direction}-regulated Genes', 
                 fontsize=13, fontweight='bold', pad=15)
    
    ax.grid(True, axis='x', alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    
    return fig

# Create individual plots
for direction in ['UP', 'DOWN']:
    fig = create_individual_barplot(direction)
    if fig is not None:
        output_file = RESULTS_FIGURES / f"panel_C_KEGG_{direction.lower()}_barplot.pdf"
        fig.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {output_file.name}")

print()

################################################################################
# Create Log File
################################################################################

from datetime import datetime

log_file = LOG_DIR / f"03_KEGG_pathway_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

with open(log_file, 'w') as f:
    f.write("KEGG PATHWAY ENRICHMENT ANALYSIS LOG (g:Profiler)\n")
    f.write(f"Execution time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Organism: {ORGANISM} (Arabidopsis thaliana)\n")
    f.write(f"Threshold: |log2FC| >= {THRESHOLD}\n")
    f.write(f"Significance threshold: {SIGNIFICANCE_THRESHOLD}\n\n")
    
    f.write("Enrichment summary:\n")
    for direction in ['UP', 'DOWN']:
        if direction in enrichment_results and len(enrichment_results[direction]) > 0:
            n_pathways = len(enrichment_results[direction])
            f.write(f"  {direction}-regulated: {n_pathways} pathways\n")
        else:
            f.write(f"  {direction}-regulated: No enrichment\n")

print(f"Log saved: {log_file.name}")

################################################################################
# Finish
################################################################################

print()
print("="*75)
print("SCRIPT 03 COMPLETED SUCCESSFULLY!")
print("="*75)
print()
print("Files created:")
print(f"  - 1 Excel file with KEGG enrichment results (results/tables/)")
print(f"  - 1 diverging barplot (recommended for manuscript)")
print(f"  - 2 individual barplots (UP and DOWN separate)")
print()
print("Next step: Run script 04_create_volcano_plot.R")
print()
