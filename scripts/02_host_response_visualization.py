# Title: Python script for DGE visualization and Pathway Analysis
# Author: Jules
# Date: 2025-08-11
#
# Description: This script reads the DGE results from the DESeq2 R script,
# generates a volcano plot, and performs GO/KEGG and GSEA pathway analysis.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
import os

print("Starting visualization and pathway analysis script...")

# --- 1. Configuration ---
DGE_RESULTS_FILE = "results/tables/DGE_LPS_vs_Control_full_results.csv"
FIGURES_DIR = "results/figures/host_response"
TABLES_DIR = "results/tables"

# Gene sets for GSEA and enrichment analysis
GENE_SETS_ENRICHMENT = ["GO_Biological_Process_2021", "KEGG_2021_Human"]
GENE_SETS_GSEA = ["KEGG_2021_Human"] # KEGG is often good for GSEA

# Significance thresholds
PADJ_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 1.0

# Create output directories
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)

# --- 2. Load DGE Results ---
print(f"Loading DGE results from {DGE_RESULTS_FILE}...")
if not os.path.exists(DGE_RESULTS_FILE):
    print(f"Error: DGE results file not found at {DGE_RESULTS_FILE}")
    print("Please run the 'scripts/run_deseq2.R' script first.")
    exit()

df = pd.read_csv(DGE_RESULTS_FILE, index_col=0)
df.dropna(inplace=True) # Remove rows with NA values, especially in padj
print("DGE results loaded successfully.")
print(df.head())

# --- 3. Generate Volcano Plot (Figure 1A) ---
print("Generating Volcano Plot...")
plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(8, 7))

# Add 'regulation' column for coloring
df['regulation'] = 'Not Significant'
df.loc[(df['padj'] < PADJ_THRESHOLD) & (df['log2FoldChange'] > LOG2FC_THRESHOLD), 'regulation'] = 'Upregulated'
df.loc[(df['padj'] < PADJ_THRESHOLD) & (df['log2FoldChange'] < -LOG2FC_THRESHOLD), 'regulation'] = 'Downregulated'

# Create scatter plot
sns.scatterplot(
    data=df,
    x='log2FoldChange',
    y=-np.log10(df['padj']),
    hue='regulation',
    palette={'Upregulated': 'red', 'Downregulated': 'blue', 'Not Significant': 'grey'},
    s=20,
    alpha=0.6,
    ax=ax,
    edgecolor=None
)

# Add threshold lines
ax.axhline(y=-np.log10(PADJ_THRESHOLD), color='k', linestyle='--', linewidth=1)
ax.axvline(x=LOG2FC_THRESHOLD, color='k', linestyle='--', linewidth=1)
ax.axvline(x=-LOG2FC_THRESHOLD, color='k', linestyle='--', linewidth=1)

# Add labels to top genes
top_genes = df[df['padj'] < PADJ_THRESHOLD].sort_values('padj').head(10)
for i, gene in top_genes.iterrows():
    ax.text(gene['log2FoldChange'], -np.log10(gene['padj']), i, fontsize=8, ha='center', va='bottom')

ax.set_title('Volcano Plot of LPS vs. Control in Microglia', fontsize=16, fontweight='bold')
ax.set_xlabel('log2(Fold Change)', fontsize=12)
ax.set_ylabel('-log10(Adjusted p-value)', fontsize=12)
ax.legend(title='Regulation', frameon=False)

volcano_path = os.path.join(FIGURES_DIR, "Fig1A_Volcano_Plot.png")
plt.savefig(volcano_path, dpi=300, bbox_inches='tight')
print(f"Volcano plot saved to {volcano_path}")
plt.close()


# --- 4. Pathway Enrichment Analysis (Figure 1B) ---
print("Running GO/KEGG Pathway Enrichment Analysis...")
# Get list of significantly upregulated genes
sig_up_genes = df[(df['regulation'] == 'Upregulated')].index.tolist()

if len(sig_up_genes) > 10:
    enr = gp.enrichr(gene_list=sig_up_genes,
                     gene_sets=GENE_SETS_ENRICHMENT,
                     organism='Human',
                     outdir=os.path.join(TABLES_DIR, 'enrichment_analysis'),
                     cutoff=0.05)

    if enr and not enr.results.empty:
        # Plotting the results (Figure 1B)
        enrichment_plot_path = os.path.join(FIGURES_DIR, "Fig1B_Enrichment_Plot.png")
        ax = gp.dotplot(enr.results,
                        title='Top Enriched Pathways in Upregulated Genes',
                        cmap='viridis_r',
                        ofname=enrichment_plot_path,
                        top_term=20)
        print(f"Enrichment plot saved to {enrichment_plot_path}")
    else:
        print("Enrichment analysis did not yield significant results.")
else:
    print("Not enough significant genes for enrichment analysis.")


# --- 5. GSEA Analysis (Figure 1C) ---
print("Running GSEA Analysis...")
# Prepare ranked list of genes for GSEA
ranked_genes = df[['log2FoldChange']].sort_values('log2FoldChange', ascending=False)
ranked_genes.reset_index(inplace=True)
ranked_genes.columns = ['gene', 'log2FoldChange']

if not ranked_genes.empty:
    pre_res = gp.prerank(rnk=ranked_genes,
                         gene_sets=GENE_SETS_GSEA,
                         organism='Human',
                         outdir=os.path.join(TABLES_DIR, 'gsea_analysis'),
                         min_size=5,
                         max_size=500,
                         permutation_num=100) # Use 100 permutations for speed in this demo

    # Plot the top result (Figure 1C)
    top_term = list(pre_res.results.keys())[0]
    gsea_plot_path = os.path.join(FIGURES_DIR, f"Fig1C_GSEA_Plot_{top_term.replace(' ', '_')}.png")
    gp.gseaplot(rank_metric=pre_res.rank_metric,
                term=top_term,
                **pre_res.results[top_term],
                ofname=gsea_plot_path)
    print(f"GSEA plot for top term '{top_term}' saved to {gsea_plot_path}")
else:
    print("Could not generate ranked gene list for GSEA.")

print("\nScript finished.")
print(f"All outputs can be found in the '{FIGURES_DIR}' and '{TABLES_DIR}' directories.")
