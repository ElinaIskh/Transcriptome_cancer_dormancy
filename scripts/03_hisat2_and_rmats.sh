#!/bin/bash

# ==============================================================================
# Project: Transcriptomic Landscape of Cancer Cell Dormancy with Proteomic Integration
# Script: 03_hisat2_and_rmats.sh
# Description: Reference genome preparation, HISAT2 genomic alignment, 
#              and Alternative Splicing analysis via rMATS.
# ==============================================================================


# --- 1. Environment Setup ---
THREADS=8
REF_DIR="./reference"
TRIMMED_DIR="./results/trimmed_data"
ALIGN_DIR="./results/alignment"
RMATS_DIR="./results/rmats_output"
TMP_DIR="./tmp"

INDEX_PREFIX="$REF_DIR/GRCh38_index"
GTF_FILE="$REF_DIR/gencode.v44.annotation.gtf"
FASTA_FILE="$REF_DIR/GRCh38.primary_assembly.genome.fa"

mkdir -p "$REF_DIR" "$ALIGN_DIR" "$RMATS_DIR" "$TMP_DIR"


# --- 2. Download and Prepare Reference (GENCODE Release 44 GRCh38) ---
if [ ! -f "$GTF_FILE" ]; then
    echo "GTF annotation missing. Downloading from GENCODE..."
    wget -c -O "${GTF_FILE}.gz" https://ebi.ac.uk
    echo "Decompressing GTF..."
    gunzip -f "${GTF_FILE}.gz"
else
    echo "GTF annotation found."
fi

if [ ! -f "$FASTA_FILE" ]; then
    echo "Genome FASTA missing. Downloading from GENCODE..."
    wget -c -O "${FASTA_FILE}.gz" https://ebi.ac.uk
    echo "Decompressing Genome FASTA (this may take a while)..."
    gunzip -f "${FASTA_FILE}.gz"
else
    echo "Genome FASTA found."
fi


# --- 3. Build HISAT2 Index ---
if [ ! -f "${INDEX_PREFIX}.1.ht2" ]; then
    echo "Step 2: Building HISAT2 index..."
    hisat2-build -p "$THREADS" "$FASTA_FILE" "$INDEX_PREFIX"
else
    echo "Step 2: HISAT2 index already exists. Skipping indexing."
fi

# --- 4. Alignment & Sorting Loop ---
echo "Step 3: Aligning reads to the genome..."

for r1 in "$TRIMMED_DIR"/*_trimmed_R1.fastq.gz; do
    [ -e "$r1" ] || continue
    
    r2="${r1/_trimmed_R1.fastq.gz/_trimmed_R2.fastq.gz}"
    base_name=$(basename "$r1" _trimmed_R1.fastq.gz)
    
    # Make name without suffixes
    sample_id=$(echo "$base_name" | awk -F'_' '{print $1"_"$2}')
    
    echo "Processing alignment for sample: $sample_id"
    
    hisat2 -p "$THREADS" \
        -x "$INDEX_PREFIX" \
        -1 "$r1" \
        -2 "$r2" \
        --summary-file "$ALIGN_DIR/${sample_id}_alignment_summary.txt" \
    | samtools sort -@ "$THREADS" -o "$ALIGN_DIR/${sample_id}.bam"
    
    samtools index "$ALIGN_DIR/${sample_id}.bam"
    samtools flagstat "$ALIGN_DIR/${sample_id}.bam" > "$ALIGN_DIR/${sample_id}_flagstat.txt"
done

# --- 5. Alternative Splicing Analysis (rMATS) ---
echo "Step 4: Preparing BAM lists and running rMATS..."

paste -sd, <(ls "$ALIGN_DIR"/*_cntrl.bam) > control_bams.txt
paste -sd, <(ls "$ALIGN_DIR"/*_exp.bam) > exp_bams.txt

echo "Control BAMs: $(cat control_bams.txt)"
echo "Experimental BAMs: $(cat exp_bams.txt)"

# Params: --readLength 150 и --variable-read-length 
rmats.py \
    --b1 control_bams.txt \
    --b2 exp_bams.txt \
    --gtf "$GTF_FILE" \
    -t paired \
    --readLength 150 \
    --variable-read-length \
    --nthread "$THREADS" \
    --od "$RMATS_DIR" \
    --tmp "$TMP_DIR"

rm -rf "$TMP_DIR"

echo "rMATS results are available in: $RMATS_DIR"
