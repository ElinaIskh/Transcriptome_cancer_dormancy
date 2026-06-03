# ==============================================================================
# Project: Integrative Transcriptomic and Proteomic Analysis of Cancer Dormancy
# Script: 06_splicing_analysis.R
# Description: Alternative splicing analysis using 'maser' on rMATS data. 
#              Includes functional profiles of splicing events and integrative 
#              comparison with DESeq2 differential expression results.
# ==============================================================================

# --- 1. Load Required Libraries ---
library(maser)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(here)

# Create output folders if they do not exist
mkdir_if_missing <- function(path) if(!dir.exists(path)) dir.create(path, recursive = TRUE)
mkdir_if_missing(here("results", "alternative_splicing"))
mkdir_if_missing(here("plots"))

# --- 2. Import rMATS data through maser ---
cat("[INFO] Importing rMATS results via maser...\n")
rmats_path <- here("results", "rmats_output")

if(!dir.exists(rmats_path)) {
  stop("Error: rMATS output directory not found! Please run the HISAT2/rMATS pipeline first.")
}

# ftype = "JC" - only Junction Counts
m <- maser(rmats_path, ftype = "JC", cond_labels = c("control", "exp"))

# Filter FDR < 0.05, |deltaPSI| > 0.1
m_sig <- topEvents(m, fdr = 0.05, deltaPSI = 0.1)
cat("[INFO] Significant splicing events identified:\n")
print(m_sig)


# --- 3. Visualisation of splicing quality (only SE type) ---
cat("[INFO] Generating diagnostic splicing plots...\n")

# Save  diagnostic splicing plots
png(here("plots", "splicing_se_volcano.png"), width = 2100, height = 1800, res = 300)
volcano(m_sig, type = "SE")
dev.off()

png(here("plots", "splicing_se_pca.png"), width = 2100, height = 1800, res = 300)
pca(m_sig, type = "SE")
dev.off()


# --- 4. Data filtration (Skipped vs Included) ---
cat("[INFO] Processing Skipped Exons (SE) details...\n")
se_results <- summary(m_sig, type = "SE")

mtor_splicing <- se_results[se_results$geneSymbol == "MTOR", ]
if(nrow(mtor_splicing) > 0) {
  cat("[INFO] Found alternative splicing events in MTOR gene!\n")
  write.table(mtor_splicing, file = here("results", "alternative_splicing", "mtor_splicing_events.tsv"), sep = "\t", row.names = FALSE)
}

# Save full table of significant SEs
write.csv(se_results, here("results", "alternative_splicing", "splicing_se_results.csv"), row.names = FALSE)

# Divide genes: skopped - (Delta PSI < -0.1) or included - (Delta PSI > 0.1)
down_genes <- unique(se_results$geneSymbol[se_results$IncLevelDifference < -0.1])
up_genes <- unique(se_results$geneSymbol[se_results$IncLevelDifference > 0.1])


# --- 5. Comparative functional analysis for splicing groups (GO and KEGG) ---
cat("[INFO] Running comparative functional analysis for splicing groups...\n")

# 5.1 Comparative GO (Biological Process)
splicing_list <- list(Skipped_Exons = down_genes, Included_Exons = up_genes)

compare_splicing_go <- compareCluster(
  geneClusters = splicing_list, 
  fun          = "enrichGO", 
  OrgDb        = org.Hs.eg.db, 
  keyType      = "SYMBOL", 
  ont          = "BP",
  pvalueCutoff = 0.05
)

if(!is.null(compare_splicing_go)) {
  splicing_go_plot <- dotplot(compare_splicing_go, showCategory = 10) + 
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 45)) +
    labs(title = "GO BP: Splicing Shifts (Delta PSI)")
  ggsave(here("plots", "splicing_go_compare.png"), plot = splicing_go_plot, width = 9, height = 8, dpi = 300)
}

# 5.2 KEGG
# COnvert Gene Symbols into Entrez ID for KEGG
ids_list <- list(
  Skipped = bitr(down_genes, "SYMBOL", "ENTREZID", "org.Hs.eg.db", drop = TRUE)$ENTREZID,
  Included = bitr(up_genes, "SYMBOL", "ENTREZID", "org.Hs.eg.db", drop = TRUE)$ENTREZID
)

compare_splicing_kegg <- compareCluster(geneClusters = ids_list, fun = "enrichKEGG", organism = "hsa", pvalueCutoff = 0.05)

if(!is.null(compare_splicing_kegg)) {
  compare_splicing_kegg_read <- setReadable(compare_splicing_kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  write.table(as.data.frame(compare_splicing_kegg_read), file = here("results", "alternative_splicing", "splicing_kegg_results.tsv"), sep = "\t", row.names = FALSE)
  
  kegg_splicing_plot <- dotplot(compare_splicing_kegg) + labs(title = "KEGG Pathways: Splicing Changes (Delta PSI)")
  ggsave(here("plots", "splicing_kegg_compare.png"), plot = kegg_splicing_plot, width = 8, height = 6, dpi = 300)
}

# 5.3 PSI visualisation for the most significant genes
if(nrow(se_results) > 0) {
  top_gene_name <- se_results$geneSymbol[which.min(se_results$PValue)]
  cat(paste("[INFO] Plotting PSI for top splicing gene:", top_gene_name, "\n"))
  
  png(here("plots", paste0("splicing_plot_gene_", top_gene_name, ".png")), width = 1800, height = 1500, res = 300)
  plotGenePSI(geneEvents(m_sig, geneS = top_gene_name), type = "SE")
  dev.off()
}


# --- 6. Integrative analysis: expression data (DESeq2) vs splicing (rMATS) ---
cat("[INFO] Integrating expression data with splicing networks...\n")

deseq_path <- here("results", "differential_expression", "differential_expression_results.tsv")

if(!file.exists(deseq_path)) {
  cat("[WARNING] 'differential_expression_results.tsv' not found. Skipping integrative integration step.\n")
} else {
  # Download DESeq2
  deseq_df <- read_delim(deseq_path, delim = "\t", show_col_types = FALSE)
  
  de_genes <- unique(deseq_df$gene_name[deseq_df$padj < 0.05 & abs(deseq_df$log2FoldChange) > 1])
  as_genes <- unique(annotation(m_sig, type = "SE")$geneSymbol)
  
  # Dual Impact Genes
  intersect_genes <- intersect(de_genes, as_genes)
  cat(paste("[INFO] Dual Impact genes identified (DEG + Splicing):", length(intersect_genes), "\n"))
  
  if(length(intersect_genes) > 0) {
    writeLines(intersect_genes, here("results", "alternative_splicing", "dual_impact_genes.txt"))
  }
  
  # Create clusters for functional analysis
  gene_clusters <- list(
    Expression  = de_genes, 
    Splicing    = as_genes, 
    Dual_Impact = intersect_genes
  )
  
  # GO-analysis(Biological Process) between expression and splicing
  compare_go_integrative <- compareCluster(
    geneClusters = gene_clusters, 
    fun          = "enrichGO", 
    OrgDb        = org.Hs.eg.db, 
    keyType      = "SYMBOL",
    ont          = "BP",
    pvalueCutoff = 0.05
  )
  
  if(!is.null(compare_go_integrative)) {
    write.table(as.data.frame(compare_go_integrative), file = here("results", "alternative_splicing", "integrative_go_expression_vs_splicing.tsv"), sep = "\t", row.names = FALSE)
    
    # Visualisation
    pdf(here("plots", "GO_genome_transcriptome.pdf"), width = 10, height = 12)
    p <- dotplot(compare_go_integrative, showCategory = 10) + 
      ggtitle("Effect of mTOR inhibitors on expression and splicing") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)
    dev.off()
    
    cat("[INFO] Integrative analysis visualization saved to 'plots/GO_genome_transcriptome.pdf'\n")
  }
}

cat("[INFO] Splicing pipeline successfully finished!\n")
