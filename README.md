# Transcriptomic Landscape of Cancer Cell Dormancy with Proteomic Integration

Repository for the final project in Bioinformatics Institute (2025-2026) and part of Master's thesis in Saint-Petersburg University (2024-2026)

**Author**: Elina Iskhakova 
> iskhakova.elina18@gmail.com
> Telegram: @Linailinai

**Supervisor**: Irina Suvorova, _Institute of Cytology of the Russian Academy of Sciences, St. Petersburg, Russia_

---
## About the project

Cancer cell dormancy is a critical latent state where tumor cells remain quiescent, often serving as a primary driver of clinical recurrence. While reversible cell cycle arrest in the G0/G1 phase is recognized as a hallmark of dormancy, the specific signaling networks and molecular markers defining this state remain poorly understood.

A well-established molecular feature of cancer dormancy across various malignancies is the suppression of the mTOR signaling pathway. Inhibiting mTOR downregulates global protein synthesis, which fundamentally reshapes the cellular proteome and rewires translational mechanisms to facilitate the dormancy transition. This proteostatic shift involves the systematic elimination of proliferation-associated proteins and the selective synthesis of proteins required for survival and adaptation.

Despite extensive research on mTOR inhibitors, global proteomic shifts and specific translational alterations in mTOR-deficient dormant cells remain largely unexplored. Furthermore, transcriptional changes alone cannot accurately predict phenotypic adaptation. To overcome this limitation, this project leverages an integrative multi-omics approach, combining RNA-Seq (transcriptomics) and Mass Spectrometry (proteomics) to provide a comprehensive, high-resolution map of the molecular landscape driving cancer cell dormancy.

## Aim and oblectives
**Aim**: To perform an integrative multi-omics analysis to elucidate the transcriptomic and proteomic landscape of cancer cells during the transition into dormancy induced by mTOR inhibition.

**Objectives**:
1. Characterize differential gene expression and key transcriptomic signatures of dormant cancer cells following mTOR inhibition.
2. Evaluate alterations in alternative splicing patterns associated with the dormancy phenotype. 
3. Profile the proteome of dormant cancer cells using mass spectrometry data to identify differentially abundant proteins. 
4. Perform integrative analysis to correlate transcriptomic shifts with proteomic data.

## Structure of repository


## Raw data
Since the raw FASTQ data is proprietary and stored locally, you need to place your raw `.fastq.bz2` files into the `raw_data/` directory before running the pipeline.

The reference transcriptome file `gencode.v44.transcripts.fa` is not included in this repository due to size constraints. You must download it directly from the [GENCODE website](https://gencodegenes.org) and place it in the project root directory before running the script.

## Workflow
### 1. RNA-seq Data Preprocessing

The preprocessing pipeline is automated via a Bash script and performs the following steps:
1. **Decompression & Conversion**: Converts `.fastq.bz2` files to `.fastq.gz` using parallel processing (`pbzip2` and `pigz`).
2. **Quality Control**: Runs `FastQC` to evaluate the initial quality of the raw sequenced reads.
3. **Trimming**: Uses `fastp` to perform automated adapter detection, sliding-window quality trimming (cutting front/tail if mean quality < 20), and filtering out short reads (< 36 bp).

### 2. Transcriptome Pseudoalignment & Quantification

To quantify transcript abundances, we use **Kallisto** for fast pseudoalignment. This step maps the trimmed paired-end reads to the reference transcriptome and generates estimated count matrices.

The script performs two main steps:
1. **Index Generation**: Builds a Kallisto index (`transcripts.idx`) from the reference FASTA file (`gencode.v44.transcripts.fa`), if it does not already exist.
2. **Abundance Quantification**: Runs `kallisto quant` in a loop across all processed samples. 

Each sample produces an independent directory containing an `abundance.tsv` file, where transcript abundances are reported in **estimated counts** and **TPM** (Transcripts Per Million). These individual tables are subsequently imported and aggregated downstream in R (using packages like `tximport`).

### Genome Alignment & Alternative Splicing Analysis

To study post-transcriptional alterations during the transition into cancer dormancy, we evaluate changes in alternative splicing. Genome alignment is performed using **HISAT2** against the **GRCh38** human genome assembly (GENCODE Release 44). Downstream alternative splicing events are quantified using **rMATS (RNA-Seq Multivariate Analysis of Transcript Splicing)**, which is the gold standard tool for this task.

### Splicing Events Evaluated
rMATS comprehensively detects and analyzes five primary types of alternative splicing events:
* **SE**: Skipped Exon
* **RI**: Retained Intron
* **MXE**: Mutually Exclusive Exons
* **A5SS**: Alternative 5' Splice Site
* **A3SS**: Alternative 3' Splice Site

The pipeline automatically handles variable read lengths (`--variable-read-length`) with a baseline reference mapping length of 150 bp for paired-end resolution.




## Usage
To execute the bash code, run the following command from the project root:
```bash
chmod +x scripts/preprocessing/SCRIPT
./scripts/preprocessing/SCRIPT
```

Before running rMATS, ensure you have the required Conda environment configured and activated:

```bash
# Create and activate rMATS environment
conda create -n rmats_env python=3.12 -y
conda activate rmats_env
conda install -c bioconda rmats -y
```

### Dependencies
Ensure the following tools are installed and available in your `PATH`:
* `pbzip2` (v1.1.13)
* `FastQC` (v0.11.9 or higher)
* `fastp` (v0.23.2 or higher v1.0.1)
* kallisto                     0.51.1
* `HISAT2` (v2.2.1 or higher)
* `Samtools` (v1.13 or higher)
* `rMATS` (v4.1.2 or higher)
