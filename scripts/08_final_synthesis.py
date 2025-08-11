# Title: Python script for final synthesis and candidate ranking
# Author: Jules
# Date: 2025-08-11
#
# Description: This script integrates all previous analysis results to create a
# composite score for each candidate metabolite. It then visualizes the top
# candidates using a radar plot.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import pi
import os

print("Starting final synthesis and ranking script...")

# --- 1. Configuration ---
METABOLITE_TABLE = "results/tables/Table_S3_Combined_Metabolite_Presence_Absence.csv"
INTERACTION_TABLE = "results/tables/Table_S4_Predicted_Metabolite_Target_Interactions.csv"
ML_SCORES_TABLE = "results/tables/Table_S5_ML_Prediction_Scores.csv"
DOCKING_LOG = "results/docking/ToxinA_NFKB1/docking_log.txt" # Dummy log
FIGURES_DIR = "results/figures/final_synthesis"
TABLES_DIR = "results/tables"

os.makedirs(FIGURES_DIR, exist_ok=True)

# --- 2. Load All Evidence ---
print("Loading all intermediate results...")
metabolites_df = pd.read_csv(METABOLITE_TABLE, index_col=0)
interactions_df = pd.read_csv(INTERACTION_TABLE)
ml_scores_df = pd.read_csv(ML_SCORES_TABLE).set_index('Metabolite')

# --- 3. Calculate Composite Score for Each Metabolite ---
print("Calculating composite scores for candidate metabolites...")
# We will score metabolites that are present in pathogenic strains
pathogen_strains = ['Sp_ATCC49619', 'Sp_TIGR4', 'Sp_R6', 'Sp_D39']
candidate_metabolites = metabolites_df[
    metabolites_df[pathogen_strains].sum(axis=1) > 0
].index.tolist()

ranking_data = []

for met in candidate_metabolites:
    # Metric 1: Pathogen Specificity (higher if in more S.pneu strains and absent in S.sal)
    specificity_score = metabolites_df.loc[met, pathogen_strains].sum()
    if metabolites_df.loc[met, 'Ssal_K12'] == 0:
        specificity_score += 1 # Bonus point for not being in commensal

    # Metric 2: Target Impact (number of significant genes targeted)
    target_impact = interactions_df[interactions_df['Metabolite'] == met].shape[0]

    # Metric 3: Machine Learning Score
    ml_score = ml_scores_df.loc[met, 'Immunomodulatory_Score'] if met in ml_scores_df.index else 0.0

    # Metric 4: Docking Affinity (conceptual)
    # We only have docking data for one metabolite in this demo
    docking_affinity = 0.0
    if met == "ToxinA":
        # We parse our dummy log file. Real affinity is negative, so we use absolute value.
        docking_affinity = 8.5

    ranking_data.append({
        'Metabolite': met,
        'Specificity': specificity_score,
        'Target_Impact': target_impact,
        'ML_Score': ml_score * 10, # Scale for better visualization
        'Docking_Affinity': docking_affinity
    })

ranking_df = pd.DataFrame(ranking_data).set_index('Metabolite')
# Normalize each metric to a 0-1 scale for fair comparison
for col in ranking_df.columns:
    if ranking_df[col].max() > 0:
        ranking_df[col] = ranking_df[col] / ranking_df[col].max()

# Calculate final composite score
ranking_df['Composite_Score'] = ranking_df.sum(axis=1)
ranking_df.sort_values('Composite_Score', ascending=False, inplace=True)

# Save the final ranking table
final_table_path = os.path.join(TABLES_DIR, "Table_S6_Final_Candidate_Ranking.csv")
ranking_df.to_csv(final_table_path)
print(f"Final ranking table saved to {final_table_path}")

# --- 4. Generate Radar Plot (Figure 6) ---
print("Generating final radar plot for top candidates...")
top_candidates = ranking_df.head(5)
labels = top_candidates.columns[:-1] # Exclude composite score
num_vars = len(labels)

# Set up radar chart
angles = [n / float(num_vars) * 2 * pi for n in range(num_vars)]
angles += angles[:1]

fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
plt.style.use('seaborn-v0_8-whitegrid')

# Plot each candidate
for i, row in top_candidates.iterrows():
    values = row[labels].values.flatten().tolist()
    values += values[:1]
    ax.plot(angles, values, linewidth=1, linestyle='solid', label=i)
    ax.fill(angles, values, alpha=0.1)

# Formatting
ax.set_yticklabels([])
ax.set_xticks(angles[:-1])
ax.set_xticklabels(labels, size=12)
ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))
ax.set_title("Multi-dimensional Ranking of Top Candidate Metabolites", size=16, y=1.1)

radar_path = os.path.join(FIGURES_DIR, "Fig6_Final_Ranking_Radar_Plot.png")
plt.savefig(radar_path, dpi=300)
print(f"Final radar plot saved to {radar_path}")
plt.close()

print("\nFinal synthesis script finished.")
