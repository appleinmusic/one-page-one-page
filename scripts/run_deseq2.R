# Title: R script for Differential Gene Expression Analysis using DESeq2
# Author: Jules
# Date: 2025-08-11
#
# Description: This script takes a raw count matrix and metadata file as input,
# performs DGE analysis using DESeq2 to compare LPS-treated vs. control samples,
# and outputs the full results table. It also creates dummy data for self-contained execution.

# --- 1. Load Libraries ---
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
if (!requireNamespace("DESeq2", quietly = TRUE))
    BiocManager::install("DESeq2")
if (!requireNamespace("readr", quietly = TRUE))
    install.packages("readr", repos = "http://cran.us.r-project.org")

library(DESeq2)
library(readr)

# --- 2. Configuration ---
counts_file <- "data/counts_matrix.tsv"
metadata_file <- "data/metadata.csv"
output_dir <- "results/tables"
output_results_file <- file.path(output_dir, "DGE_LPS_vs_Control_full_results.csv")

# Create directories if they don't exist
dir.create("data", showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- 3. Create Dummy Data for Demonstration ---
# In a real scenario, these files would be the output of an RNA-seq quantification pipeline (e.g., salmon, STAR)
# Dummy counts matrix (1000 genes, 6 samples)
set.seed(42)
counts_matrix <- data.frame(
  ctrl1 = rnbinom(1000, mu = 500, size = 10),
  ctrl2 = rnbinom(1000, mu = 500, size = 10),
  ctrl3 = rnbinom(1000, mu = 500, size = 10),
  lps1 = rnbinom(1000, mu = 1500, size = 10),
  lps2 = rnbinom(1000, mu = 1500, size = 10),
  lps3 = rnbinom(1000, mu = 1500, size = 10)
)
rownames(counts_matrix) <- paste0("Gene", 1:1000)
# Introduce some actual differential expression for ~100 genes
de_genes <- sample(1:1000, 100)
counts_matrix[de_genes, 4:6] <- rnbinom(100, mu = 4000, size = 10)
write.table(counts_matrix, counts_file, sep="\t", quote=FALSE, col.names=NA)

# Dummy metadata
metadata <- data.frame(
  SampleID = c("ctrl1", "ctrl2", "ctrl3", "lps1", "lps2", "lps3"),
  Condition = c("Control", "Control", "Control", "LPS", "LPS", "LPS")
)
write.csv(metadata, metadata_file, quote=FALSE, row.names=FALSE)


# --- 4. Load Data ---
# Load count data
count_data_df <- read.table(counts_file, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
# Load metadata
col_data_df <- read.csv(metadata_file, row.names=1)

# Ensure row names of metadata match column names of count data
if (!all(rownames(col_data_df) %in% colnames(count_data_df))) {
    stop("Mismatch between metadata rows and count data columns. Please check sample IDs.")
}
# Order count matrix columns to match metadata rows
count_data_ordered <- count_data_df[, rownames(col_data_df)]


# --- 5. Create DESeqDataSet Object ---
# The design formula `~ Condition` tells DESeq2 to model the counts based on the 'Condition' column
dds <- DESeqDataSetFromMatrix(countData = count_data_ordered,
                              colData = col_data_df,
                              design = ~ Condition)

# --- 6. Pre-filtering ---
# Keep only rows that have at least 10 reads total across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# --- 7. Run DESeq Analysis ---
# Set the reference level (the group to compare against)
dds$Condition <- relevel(dds$Condition, ref = "Control")

# Run the DESeq function
dds <- DESeq(dds)

# --- 8. Get Results ---
# Get the results for the comparison of 'LPS' vs 'Control'
res <- results(dds, name="Condition_LPS_vs_Control")

# Sort results by adjusted p-value
res_ordered <- res[order(res$padj),]

# --- 9. Save Results ---
# Convert to a data frame and save
write.csv(as.data.frame(res_ordered), file = output_results_file)

cat("DESeq2 analysis complete. Results saved to:", output_results_file, "\n")
