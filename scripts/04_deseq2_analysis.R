# ==============================================================================
# Project: Transcriptomic Landscape of Cancer Cell Dormancy with Proteomic Integration
# Script: 04_deseq2_analysis.R
# Description: Differential Gene Expression (DGE) analysis with DESeq2. 
#              Includes multi-factor design (controlling for cell line batch effect),
#              quality control (PCA, sample distance heatmap), and annotation.
# ==============================================================================


# 1. Load Required Libraries
library(rtracklayer)
library(tximport)
library(DESeq2)
library(rhdf5)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(here) # For platform-independent, relative paths

# Create output folders if they do not exist
mkdir_if_missing <- function(path) if(!dir.exists(path)) dir.create(path, recursive = TRUE)
mkdir_if_missing(here("results", "tables"))
mkdir_if_missing(here("plots"))


# 2. Experimental Design Setup
samples <- c("A549_cntrl", "A549_exp", "T98G_cntrl", "T98G_exp", "PA1_cntrl", "PA1_exp")
condition <- factor(c("control", "treated", "control", "treated", "control", "treated"))
cell_line <- factor(c("A549", "A549", "T98G", "T98G", "PA1", "PA1"))

colData_cell <- data.frame(
  row.names = samples,
  cell_line = cell_line,
  condition = condition
)

# 3. Dynamic Path Construction for Kallisto Abundance Files
# Looks for results/kallisto_output/[sample_name]/abundance.h5
files <- here("results", "kallisto_output", samples, "abundance.h5")
names(files) <- samples

if(!all(file.exists(files))) {
  stop("Error: One or more abundance.h5 files are missing in results/kallisto_output/!")
}


# 4. Import Transcriptome Annotation & Build tx2gene Map
gtf_path <- here("reference", "gencode.v44.annotation.gtf")
if(!file.exists(gtf_path)) {
  stop("Error: Annotation file 'gencode.v44.annotation.gtf' not found in 'reference/' directory!")
}

echo <- function(msg) cat(paste0("[INFO] ", msg, "\n"))
echo("Importing GTF annotation...")
gtf <- import(gtf_path)

tx2gene <- unique(data.frame(
  transcript_id = gtf$transcript_id,
  gene_id = gtf$gene_id
))


# 5. Data Import via tximport
echo("Importing Kallisto quantification data...")
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)

# Quick sanity check on imported dimensions
echo(paste("Imported counts matrix shape:", paste(dim(txi$counts), collapse = "x")))


# 6. DESeq2 Analysis (Controlling for Cell Line Batch Effect)
echo("Initializing DESeq2 analysis...")
dds_cell <- DESeqDataSetFromTximport(
  txi,
  colData = colData_cell,
  design = ~ cell_line + condition # Multi-factor design to control for biological variance between cell lines
)

# Set control as reference group
dds_cell$condition <- relevel(dds_cell$condition, ref = "control")

# Run differential expression pipeline
dds_cell <- DESeq(dds_cell)


# 7. Quality Control & Data Visualization
echo("Performing Variance Stabilizing Transformation (VST)...")
vsd_cell <- vst(dds_cell, blind = FALSE)

# 7.1 PCA Plot Generation
echo("Generating PCA Plot...")
pcaData <- plotPCA(vsd_cell, intgroup = c("condition", "cell_line"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = condition, shape = cell_line)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("control" = "#7570b3", "treated" = "#d95f02")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  coord_fixed() +
  labs(title = "Principal Component Analysis (PCA)")

# Save PCA plot to file
ggsave(here("plots", "pca_plot.png"), plot = pca_plot, width = 7, height = 6, dpi = 300)


# 7.2 Sample Distance Heatmap
echo("Generating Sample Distance Heatmap...")
sampleDists <- dist(t(assay(vsd_cell)))
sampleDistMatrix <- as.matrix(sampleDists)

png(here("plots", "sample_distance_heatmap.png"), width = 2100, height = 1800, res = 300)
pheatmap(
  sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  main = "Sample-to-Sample Distances"
)
dev.off()


# 8. Results Extraction and Gene Annotation
echo("Extracting and annotating results...")
res <- results(dds_cell)

# Filter out rows with NA padj values and sort by statistical significance
res <- res[!is.na(res$padj), ]
sorted_res <- res[order(res$padj), ]

sorted_df <- data.frame("gene_id" = rownames(sorted_res), sorted_res)

# Build a fast mapping dataframe for human-readable gene symbols
gene_annot <- unique(data.frame(
  gene_id = gtf$gene_id,
  gene_name = gtf$gene_name
))

# Merge tables to add 'gene_name' column
res_annot <- merge(sorted_df, gene_annot, by = "gene_id", all.x = TRUE)
new_order <- c("gene_name", "gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
res_annot <- res_annot[, new_order]

# Save final annotated table
output_table_path <- here("results", "tables", "differential_expression_results.tsv")
write.table(res_annot, file = output_table_path, sep = "\t", row.names = FALSE, quote = FALSE)

echo(paste("Analysis complete! Annotated results saved to:", output_table_path))
