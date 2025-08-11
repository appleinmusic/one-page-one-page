# Title: Python script for comparative analysis of bacterial metabolites
# Author: Jules
# Date: 2025-08-11
#
# Description: This script reads the reconstructed metabolite lists for multiple
# strains, performs a comparative analysis, and generates an Upset plot and
# a pathway heatmap to visualize the core and accessory metabolomes.

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from upsetplot import from_contents, UpSet
import os

print("Starting bacterial metabolism comparison script...")

# --- 1. Configuration ---
RECONSTRUCTION_DIR = "results/metabolic_models"
FIGURES_DIR = "results/figures/bacterial_metabolism"
TABLES_DIR = "results/tables"

STRAIN_PREFIXES = [
    "Sp_ATCC49619", "Sp_TIGR4", "Sp_R6", "Sp_D39", "Ssal_K12"
]

os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)

# --- 2. Load and Parse Metabolite Data ---
print("Loading metabolite data for all strains...")
metabolite_sets = {}
all_metabolites = set()

for prefix in STRAIN_PREFIXES:
    file_path = os.path.join(RECONSTRUCTION_DIR, f"{prefix}-all-Metabolites.tbl")
    if os.path.exists(file_path):
        # Skip header, read the first column (ID)
        df = pd.read_csv(file_path, sep='\t')
        # Use a generic metabolite ID for comparison, e.g., from 'Name' or 'ID' column
        # Here we use 'Name' for better readability
        strain_metabolites = set(df['Name'])
        metabolite_sets[prefix] = strain_metabolites
        all_metabolites.update(strain_metabolites)
        print(f"  > Loaded {len(strain_metabolites)} metabolites for {prefix}")
    else:
        print(f"Warning: Metabolite file not found for {prefix}")
        metabolite_sets[prefix] = set()

# --- 3. Generate Upset Plot (Figure 2A) ---
print("Generating Upset plot...")
# The UpSet plot requires data in a specific format.
# `from_contents` is a helper function to create this from a dictionary of sets.
upset_data = from_contents(metabolite_sets)

plt.style.use('seaborn-v0_8-whitegrid')
fig = plt.figure(figsize=(12, 7))
upset = UpSet(upset_data, subset_size='count', show_counts=True, sort_by='degree')
upset.plot(fig=fig)
plt.suptitle("Core and Accessory Metabolome of Streptococcus Strains", fontsize=16)

upset_plot_path = os.path.join(FIGURES_DIR, "Fig2A_Metabolite_Upset_Plot.png")
plt.savefig(upset_plot_path, dpi=300, bbox_inches='tight')
print(f"Upset plot saved to {upset_plot_path}")
plt.close()


# --- 4. Generate Metabolic Pathway Heatmap (Figure 2B) ---
# This is a conceptual demonstration. A real implementation would involve mapping
# each metabolite to its KEGG pathway and then checking for pathway completeness.
print("Generating conceptual pathway presence/absence heatmap...")

# Dummy pathway data: map some metabolites to pathways
pathway_mapping = {
    'D-Glucose': 'Glycolysis',
    'Pyruvate': 'Glycolysis',
    'ToxinA': 'Virulence Factors',
    'Salivaricin': 'Bacteriocin Production'
}
all_pathways = sorted(list(set(pathway_mapping.values())))
pathway_presence = pd.DataFrame(0, index=all_pathways, columns=STRAIN_PREFIXES)

for strain, metabolites in metabolite_sets.items():
    for met in metabolites:
        if met in pathway_mapping:
            pathway = pathway_mapping[met]
            pathway_presence.loc[pathway, strain] = 1

fig, ax = plt.subplots(figsize=(8, 5))
sns.heatmap(
    pathway_presence,
    annot=True,
    cmap="YlGnBu",
    linewidths=.5,
    ax=ax,
    cbar=False
)
ax.set_title("Presence of Key Metabolic Pathways Across Strains", fontsize=14, fontweight='bold')
ax.set_xlabel("Strain", fontsize=12)
ax.set_ylabel("Metabolic Pathway", fontsize=12)
plt.xticks(rotation=45, ha='right')

heatmap_path = os.path.join(FIGURES_DIR, "Fig2B_Pathway_Heatmap.png")
plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
print(f"Pathway heatmap saved to {heatmap_path}")
plt.close()


# --- 5. Save Combined Metabolite Table (Table S3) ---
print("Generating and saving combined metabolite table (Table S3)...")
# Create a presence/absence dataframe
presence_absence_df = pd.DataFrame(0, index=sorted(list(all_metabolites)), columns=STRAIN_PREFIXES)
for strain, metabolites in metabolite_sets.items():
    presence_absence_df.loc[list(metabolites), strain] = 1

table_path = os.path.join(TABLES_DIR, "Table_S3_Combined_Metabolite_Presence_Absence.csv")
presence_absence_df.to_csv(table_path)
print(f"Combined metabolite table saved to {table_path}")

print("\nScript finished.")
print(f"All outputs can be found in the '{FIGURES_DIR}' and '{TABLES_DIR}' directories.")
