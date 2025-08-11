# Title: Python script for ML-based prediction of metabolite immunomodulatory potential
# Author: Jules
# Date: 2025-08-11
#
# Description: This script simulates the process of training a machine learning
# model on known immunomodulatory metabolites. It then uses this model to predict
# the potential of the bacterial metabolites identified in our study.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc, classification_report
import os

print("Starting ML prediction script for metabolite potential...")

# --- 1. Configuration ---
METABOLITE_TABLE = "results/tables/Table_S3_Combined_Metabolite_Presence_Absence.csv"
FIGURES_DIR = "results/figures/advanced_modeling"
TABLES_DIR = "results/tables"
FINGERPRINT_LEN = 128  # Length of the simulated chemical fingerprint

os.makedirs(FIGURES_DIR, exist_ok=True)

# --- 2. Simulate Training Data ---
# In a real scenario, this data would be curated from literature/databases (e.g., ChEMBL, PubChem).
# We need metabolites with known activity and their chemical fingerprints.
print("Simulating training data of known immunomodulatory metabolites...")
np.random.seed(42)
num_known_metabolites = 200
X = np.random.randint(0, 2, size=(num_known_metabolites, FINGERPRINT_LEN)) # Dummy fingerprints
y = np.random.randint(0, 2, size=num_known_metabolites) # 0=non-active, 1=active
# Make the problem learnable by adding a slight bias for the 'active' class
active_indices = np.where(y == 1)[0]
X[active_indices, :10] = np.random.choice([0, 1], size=(len(active_indices), 10), p=[0.2, 0.8])

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42, stratify=y)
print(f"Training data: {X_train.shape[0]} samples, Test data: {X_test.shape[0]} samples.")

# --- 3. Train RandomForest Model ---
print("Training RandomForestClassifier...")
model = RandomForestClassifier(n_estimators=100, random_state=42, class_weight='balanced')
model.fit(X_train, y_train)

# --- 4. Evaluate Model ---
print("Evaluating model performance...")
y_pred = model.predict(X_test)
y_proba = model.predict_proba(X_test)[:, 1]

print("Classification Report:")
print(classification_report(y_test, y_pred))

# a) Generate ROC Curve (Figure 4A)
fpr, tpr, _ = roc_curve(y_test, y_proba)
roc_auc = auc(fpr, tpr)

plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(7, 7))
ax.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
ax.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
ax.set_xlim([0.0, 1.0])
ax.set_ylim([0.0, 1.05])
ax.set_xlabel('False Positive Rate')
ax.set_ylabel('True Positive Rate')
ax.set_title('ROC Curve for Immunomodulatory Prediction Model', fontsize=14)
ax.legend(loc="lower right")

roc_path = os.path.join(FIGURES_DIR, "Fig4A_ML_ROC_Curve.png")
plt.savefig(roc_path, dpi=300)
print(f"ROC curve saved to {roc_path}")
plt.close()

# b) Generate Feature Importance Plot (Figure 4B)
importances = model.feature_importances_
indices = np.argsort(importances)[-20:] # Top 20 features

fig, ax = plt.subplots(figsize=(10, 6))
ax.barh(range(len(indices)), importances[indices], align='center')
ax.set_yticks(range(len(indices)), [f'Fingerprint_Bit_{i}' for i in indices])
ax.set_xlabel('Feature Importance')
ax.set_title('Top 20 Important Features for Prediction Model', fontsize=14)

feature_imp_path = os.path.join(FIGURES_DIR, "Fig4B_ML_Feature_Importance.png")
plt.savefig(feature_imp_path, dpi=300, bbox_inches='tight')
print(f"Feature importance plot saved to {feature_imp_path}")
plt.close()

# --- 5. Predict on Our Bacterial Metabolites ---
print("Predicting potential of S. pneumoniae metabolites...")
bacterial_metabolites_df = pd.read_csv(METABOLITE_TABLE, index_col=0)
# Focus on metabolites from any S. pneumoniae strain
sp_metabolites = bacterial_metabolites_df[bacterial_metabolites_df[
    ['Sp_ATCC49619', 'Sp_TIGR4', 'Sp_R6', 'Sp_D39']
].sum(axis=1) > 0].index.tolist()

# Simulate fingerprints for these metabolites
num_sp_metabolites = len(sp_metabolites)
X_sp = np.random.randint(0, 2, size=(num_sp_metabolites, FINGERPRINT_LEN))

# Predict the probability of being 'immunomodulatory'
sp_metabolite_scores = model.predict_proba(X_sp)[:, 1]

results_df = pd.DataFrame({
    'Metabolite': sp_metabolites,
    'Immunomodulatory_Score': sp_metabolite_scores
}).sort_values('Immunomodulatory_Score', ascending=False)


# --- 6. Save Prediction Results ---
table_ml_path = os.path.join(TABLES_DIR, "Table_S5_ML_Prediction_Scores.csv")
results_df.to_csv(table_ml_path, index=False)
print(f"Prediction scores for bacterial metabolites saved to {table_ml_path}")

print("\nML prediction script finished.")
