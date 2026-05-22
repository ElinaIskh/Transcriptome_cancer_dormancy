#!/bin/bash

# ==============================================================================
# Project: Transcriptomic Landscape of Cancer Cell Dormancy with Proteomic Integration
# Script: 02_kallisto.sh
# Description: Transcriptome index generation and pseudoalignment/quantification
#              using Kallisto for paired-end RNA-Seq data.
# ==============================================================================


# --- 1. Environment Setup ---
# Define directory paths (adjust these according to your local environment)
REFERENCE_FA="gencode.v44.transcripts.fa"
INDEX_FILE="transcripts.idx"
TRIMMED_DATA_DIR="./results/trimmed_data"
KALLISTO_OUT_DIR="./results/kallisto_output"

# Create output directory for quantification results
mkdir -p "$KALLISTO_OUT_DIR"


# --- 2. Index Generation ---
# Generate the transcriptome index if it does not already exist
if [ ! -f "$INDEX_FILE" ]; then
    echo "Step 1: Index file '$INDEX_FILE' not found. Building index from $REFERENCE_FA..."
    if [ ! -f "$REFERENCE_FA" ]; then
        echo "Error: Reference file '$REFERENCE_FA' is missing! Please download it first."
        exit 1
    fi
    kallisto index -i "$INDEX_FILE" "$REFERENCE_FA"
    echo "Index successfully built."
else
    echo "Step 1: Index file '$INDEX_FILE' already exists. Skipping index generation."
fi


# --- 3. Quantification (kallisto quant) ---
echo "Step 2: Quantifying expression levels for all samples..."

# Loop through all forward trimmed reads and process paired-end samples
for r1 in "$TRIMMED_DATA_DIR"/*_trimmed_R1.fastq.gz; do
    # Check if files exist to avoid loop errors
    [ -e "$r1" ] || continue
    
    # Identify corresponding reverse trimmed read
    r2="${r1/_trimmed_R1.fastq.gz/_trimmed_R2.fastq.gz}"
    
    # Extract baseline sample name for creating unique output folders (e.g., S_1 or A549_cntrl)
    base_name=$(basename "$r1" _trimmed_R1.fastq.gz)
    
    # Define a sample-specific output folder
    sample_out_dir="$KALLISTO_OUT_DIR/${base_name}_kallisto"
    
    echo "Processing sample: $base_name"
    
    kallisto quant \
        -i "$INDEX_FILE" \
        -o "$sample_out_dir" \
        "$r1" \
        "$r2"
done

echo "=== Kallisto processing completed successfully! ==="
echo "Output 'abundance.tsv' files can be found in: $KALLISTO_OUT_DIR"
