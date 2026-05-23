#!/bin/bash


# ==============================================================================
# Project: Integrative Transcriptomic and Proteomic Analysis of Cancer Dormancy
# Script: 01_rna_seq_preprocessing.sh
# Description: Automated pipeline for RNA-Seq raw data preprocessing.
#              Includes decompression (.bz2 to .gz), Quality Control (FastQC), 
#              and adapter/quality trimming (fastp).
# ==============================================================================

# --- 1. Environment Setup ---
# Define directory paths (adjust these according to your local environment)
RAW_DATA_DIR="./raw_data"
FASTQC_RAW_DIR="./results/fastqc_raw"
TRIMMED_DATA_DIR="./results/trimmed_data"
FASTP_REPORTS_DIR="./results/fastp_reports"

# Create output directories if they do not exist
mkdir -p "$FASTQC_RAW_DIR" "$TRIMMED_DATA_DIR" "$FASTP_REPORTS_DIR"


# --- 2. Decompression and Format Conversion (.bz2 to .gz) ---
# Multi-threaded conversion from bzip2 to gzip for improved compatibility
echo "Step 1: Converting files from .fastq.bz2 to .fastq.gz..."
cd "$RAW_DATA_DIR"

for f in *.fastq.bz2; do
    # Check if files exist to avoid loop errors
    [ -e "$f" ] || continue
    
    output_gz="${f%.bz2}.gz"
    if [ ! -f "$output_gz" ]; then
        echo "Processing: $f -> $output_gz"
        pbzip2 -dc "$f" | pigz > "$output_gz"
    else
        echo "Skipping: $output_gz already exists."
    fi
done

cd .. # Return to the project root directory


# --- 3. Initial Quality Control (FastQC) ---
echo "Step 2: Running FastQC on raw data..."
# Run FastQC in parallel using multiple files at once
fastqc -o "$FASTQC_RAW_DIR" "$RAW_DATA_DIR"/*.fastq.gz


# --- 4. Quality and Adapter Trimming (fastp) ---
echo "Step 3: Trimming adapters and low-quality reads with fastp..."

# Loop through forward reads (_1.fastq.gz) and automatically find reverse reads (_2.fastq.gz)
for r1 in "$RAW_DATA_DIR"/*_1.fastq.gz; do
    [ -e "$r1" ] || continue
    
    # Identify corresponding reverse read
    r2="${r1/_1.fastq.gz/_2.fastq.gz}"
    
    # Extract sample baseline name for naming outputs (e.g., S_1)
    base_name=$(basename "$r1" _1.fastq.gz)
    
    echo "Trimming sample: $base_name"
    
    fastp \
        -i "$r1" \
        -I "$r2" \
        -o "$TRIMMED_DATA_DIR/${base_name}_trimmed_R1.fastq.gz" \
        -O "$TRIMMED_DATA_DIR/${base_name}_trimmed_R2.fastq.gz" \
        --detect_adapter_for_pe \
        --cut_front \
        --cut_tail \
        --cut_window_size 4 \
        --cut_mean_quality 20 \
        --length_required 36 \
        --html "$FASTP_REPORTS_DIR/${base_name}_fastp.html" \
        --json "$FASTP_REPORTS_DIR/${base_name}_fastp.json"
done
