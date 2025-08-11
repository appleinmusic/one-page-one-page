# Colab-Friendly Pipeline: Predicting S. pneumoniae Metabolite Effects on Microglia

## 1. Project Overview

This project provides a high-end, reproducible bioinformatics pipeline designed to run entirely within **Google Colab**. It addresses the research question of how *Streptococcus pneumoniae* metabolites might influence microglial activation, using a modern, API-driven approach that requires no local installation of heavy bioinformatics software.

The entire analysis is broken down into a series of Jupyter Notebooks. Each notebook is self-contained, installing its own dependencies and fetching real data directly from public scientific databases like NCBI GEO, KEGG, and STITCH.

---

## 2. How to Run This Project

Click the "Open in Colab" badges below to launch each notebook directly in your browser. **It is recommended to run them in sequential order (01 through 04).**

| #  | Notebook                               | Description                                                                 | Link                                                                                                                                                |
|----|----------------------------------------|-----------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------|
| 01 | `01_Host_Response_Analysis.ipynb`      | Defines the microglial activation signature from real RNA-seq data.         | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YOUR_GITHUB_USERNAME/YOUR_REPO_NAME/blob/main/01_Host_Response_Analysis.ipynb)      |
| 02 | `02_Bacterial_Metabolite_Analysis.ipynb` | Fetches bacterial metabolite data from the KEGG database.                   | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YOUR_GITHUB_USERNAME/YOUR_REPO_NAME/blob/main/02_Bacterial_Metabolite_Analysis.ipynb) |
| 03 | `03_Bridging_and_ML_Analysis.ipynb`      | Connects metabolites to host genes (STITCH API) & predicts their potential (ML). | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YOUR_GITHUB_USERNAME/YOUR_REPO_NAME/blob/main/03_Bridging_and_ML_Analysis.ipynb)      |
| 04 | `04_Final_Synthesis.ipynb`             | Integrates all results and ranks the top candidate metabolites.             | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YOUR_GITHUB_USERNAME/YOUR_REPO_NAME/blob/main/04_Final_Synthesis.ipynb)      |

**Note:** You will need to replace `YOUR_GITHUB_USERNAME/YOUR_REPO_NAME` with your actual GitHub repository path once you have uploaded these files.

---

## 3. The Analysis Workflow

### Notebook 01: Host Response Analysis
- **Input**: GEO Accession ID (`GSE155408`).
- **Process**: Downloads data using `GEOparse`, performs DGE analysis with `pydeseq2`, and runs pathway analysis with `gseapy`.
- **Output**: A list of significantly changed host genes and pathways that define "microglial activation".

### Notebook 02: Bacterial Metabolite Analysis
- **Input**: KEGG organism codes (e.g., `spn` for *S. pneumoniae* TIGR4).
- **Process**: Queries the KEGG REST API to get all metabolites for each selected bacterium.
- **Output**: A comparative table of metabolite presence/absence across strains.

### Notebook 03: Bridging & Machine Learning
- **Input**: Results from Notebooks 01 & 02.
- **Process**:
    1.  Queries the STITCH API to find predicted protein targets for key bacterial metabolites.
    2.  Filters these targets against the significant host gene list.
    3.  Loads a real-world dataset of immunomodulatory molecules to train a Random Forest classifier.
    4.  (Conceptually) Predicts the potential of the bacterial metabolites.
- **Output**: A network of metabolite-gene interactions and a table of ML-based scores.

### Notebook 04: Final Synthesis
- **Input**: Results from all previous notebooks.
- **Process**: Calculates a composite score for each metabolite based on multiple lines of evidence (pathogen specificity, target impact, ML score).
- **Output**: A final ranked list of candidate metabolites, visualized as a radar plot.

---
This project structure ensures full reproducibility and accessibility for users without access to dedicated bioinformatics servers.
