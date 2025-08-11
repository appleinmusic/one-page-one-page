#!/bin/bash
# Title: Shell script for bacterial genome annotation and metabolic reconstruction
# Author: Jules
# Date: 2025-08-11
#
# Description: This script orchestrates the download of bacterial genomes,
# annotation with Prokka, and metabolic network reconstruction with gapseq.
# It then calls a Python script to perform the comparative analysis.

set -e # Exit immediately if a command exits with a non-zero status.

echo "--- Stage 3: Bacterial Metabolism Pipeline ---"

# --- 1. Configuration ---
GENOME_DIR="data/genomes"
ANNOTATION_DIR="results/annotations"
RECONSTRUCTION_DIR="results/metabolic_models"

# Note: In a real run, we would use the accessions gathered in Stage 1.
# For this script, we will use placeholder names.
declare -A ACCESSIONS
ACCESSIONS=(
    ["Sp_ATCC49619"]="AP018938.1"
    ["Sp_TIGR4"]="GCA_000006885.1"
    ["Sp_R6"]="GCA_000007045.1"
    ["Sp_D39"]="GCA_000014365.2"
    ["Ssal_K12"]="GCA_049310985.1"
)

mkdir -p $GENOME_DIR $ANNOTATION_DIR $RECONSTRUCTION_DIR

# --- 2. Download Genomes (Conceptual Step) ---
echo "[1/4] Downloading genomes from NCBI..."
for prefix in "${!ACCESSIONS[@]}"; do
    acc=${ACCESSIONS[$prefix]}
    echo "  > Downloading ${prefix} (${acc})..."
    # The real command would be:
    # ncbi-datasets download genome accession $acc --filename "${GENOME_DIR}/${prefix}.zip"
    # unzip -o "${GENOME_DIR}/${prefix}.zip" -d "${GENOME_DIR}/${prefix}_files"
    # mv "${GENOME_DIR}/${prefix}_files/ncbi_dataset/data/${acc}/*.fna" "${GENOME_DIR}/${prefix}.fna"

    # For demonstration, we'll just create dummy files.
    echo "> ${prefix}" > "${GENOME_DIR}/${prefix}.fna"
done
echo "Genome download conceptually complete."
echo

# --- 3. Annotate Genomes with Prokka (Conceptual Step) ---
echo "[2/4] Annotating genomes with Prokka..."
for prefix in "${!ACCESSIONS[@]}"; do
    echo "  > Annotating ${prefix}..."
    # The real command would be:
    # prokka --outdir "${ANNOTATION_DIR}/${prefix}" --prefix "${prefix}" --locustag "${prefix}" --cpus 0 "${GENOME_DIR}/${prefix}.fna"

    # For demonstration, create a dummy output directory and the key .gbk file.
    mkdir -p "${ANNOTATION_DIR}/${prefix}"
    echo "LOCUS       ${prefix}    2000000 bp    DNA     circular BCT 11-AUG-2025" > "${ANNOTATION_DIR}/${prefix}/${prefix}.gbk"
done
echo "Annotation conceptually complete."
echo

# --- 4. Reconstruct Metabolic Models with gapseq (Conceptual Step) ---
echo "[3/4] Reconstructing metabolic models with gapseq..."
for prefix in "${!ACCESSIONS[@]}"; do
    echo "  > Running gapseq for ${prefix}..."
    # The real commands would be:
    # gapseq find -p "${prefix}" -t Bacteria -b 100 "${ANNOTATION_DIR}/${prefix}/${prefix}.gbk"
    # gapseq fill -m "${prefix}-model.RDS" -c "${prefix}-blast-out.txt" -b 100
    # gapseq extract -m "${prefix}-model.RDS" -f "${RECONSTRUCTION_DIR}"

    # For demonstration, create a dummy metabolite list file.
    # Gapseq typically produces a file like 'prefix-all-Metabolites.tbl'
    METABOLITE_FILE="${RECONSTRUCTION_DIR}/${prefix}-all-Metabolites.tbl"
    echo -e "ID\tName\tFormula\tCharge" > "${METABOLITE_FILE}"
    echo -e "cpd00001\tH2O\tH2O\t0" >> "${METABOLITE_FILE}"
    echo -e "cpd00027\tD-Glucose\tC6H12O6\t0" >> "${METABOLITE_FILE}"
    echo -e "cpd00020\tPyruvate\tC3H3O3\t-1" >> "${METABOLITE_FILE}"
    # Add some strain-specific dummy metabolites for interesting plots
    if [[ $prefix == "Sp_TIGR4" || $prefix == "Sp_D39" ]]; then
      echo -e "cpd_virulence_A\tToxinA\tC10H12N2O\t0" >> "${METABOLITE_FILE}"
    fi
    if [[ $prefix == "Ssal_K12" ]]; then
      echo -e "cpd_probiotic_B\tSalivaricin\tC20H30N5O5\t1" >> "${METABOLITE_FILE}"
    fi
done
echo "Metabolic reconstruction conceptually complete."
echo

# --- 5. Run Python script for comparative analysis and visualization ---
echo "[4/4] Running Python script for comparative analysis..."
# This assumes the python script `04_bacterial_metabolism_comparison.py` exists
# and that python is in the environment.
# We will create this script in the next step.
# python3 scripts/04_bacterial_metabolism_comparison.py

echo "--- Stage 3 Pipeline Finished (Conceptual Execution) ---"
echo "The next step is to create the Python script to process these results."
