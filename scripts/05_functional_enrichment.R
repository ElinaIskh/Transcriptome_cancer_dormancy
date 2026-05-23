# ==============================================================================
# Project: Integrative Transcriptomic and Proteomic Analysis of Cancer Dormancy
# Script: 05_functional_enrichment.R
# Description: Functional enrichment analysis (GO and KEGG) for differentially 
#              expressed genes (DEGs) using clusterProfiler.
# ==============================================================================

# 1. Load Required Libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(here)

# Create output folders if they do not exist
if(!dir.exists(here("results", "differential_expression"))) dir.create(here("results", "differential_expression"), recursive = TRUE)
if(!dir.exists(here("plots"))) dir.create(here("plots"), recursive = TRUE)

cat("[INFO] Loading differential expression results...\n")

# 2. Download differential expression results
deseq_results_path <- here("results", "differential_expression", "differential_expression_results.tsv")

if(!file.exists(deseq_results_path)) {
  stop("Error: 'differential_expression_results.tsv' not found! Please run the DESeq2 script first.")
}

res_df <- read.table(deseq_results_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 3. Filter differential expressed genes
# padj < 0.05 and |log2FC| > 1
p_threshold <- 0.05
lfc_threshold <- 1

up_genes_df <- res_df %>% filter(padj < p_threshold & log2FoldChange > lfc_threshold)
down_genes_df <- res_df %>% filter(padj < p_threshold & log2FoldChange < -lfc_threshold)
all_de_genes_df <- res_df %>% filter(padj < p_threshold & abs(log2FoldChange) > lfc_threshold)

cat(paste("[INFO] Detected DEGs - Total:", nrow(all_de_genes_df), 
          "| Up-regulated:", nrow(up_genes_df), 
          "| Down-regulated:", nrow(down_genes_df), "\n"))

# Convert ENSEMBL ID 
get_clean_ensembl <- function(gene_ids) {
  gsub("\\..*$", "", gene_ids)
}

de_genes_ensembl <- get_clean_ensembl(all_de_genes_df$gene_id)
universe_genes_ensembl <- get_clean_ensembl(res_df$gene_id)


# ==============================================================================
# 4. Gene Ontology - GO Enrichment Analysis
# ==============================================================================
cat("[INFO] Running Gene Ontology (GO) enrichment analysis...\n")

ego_bp <- enrichGO(
  gene          = de_genes_ensembl,
  universe      = universe_genes_ensembl,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP", # BP = Biological Process 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE  # Convert to Gene Symbols
)

# Save GO results as table
if(!is.null(ego_bp) && nrow(ego_bp) > 0) {
  write.table(as.data.frame(ego_bp), 
              file = here("results", "differential_expression", "go_bp_enrichment_results.tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Visualise GO: Dotplot
  go_dotplot <- dotplot(ego_bp, showCategory = 20) + 
    labs(title = "GO Biological Process Enrichment")
  ggsave(here("plots", "go_bp_dotplot.png"), plot = go_dotplot, width = 9, height = 8, dpi = 300)
  
  # Visualise GO: cnetplot
  go_cnet <- cnetplot(ego_bp, showCategory = 5, foldChange = NULL) + 
    labs(title = "GO Gene-Concept Network")
  ggsave(here("plots", "go_bp_cnetplot.png"), plot = go_cnet, width = 10, height = 8, dpi = 300)
} else {
  cat("[WARNING] No significant GO terms found.\n")
}


# ==============================================================================
# 5. Signaling pathways analysis: KEGG
# ==============================================================================
cat("[INFO] Running KEGG pathway enrichment analysis...\n")

gene_conversion <- bitr(
  de_genes_ensembl, 
  fromType = "ENSEMBL", 
  toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db
)

universe_conversion <- bitr(
  universe_genes_ensembl, 
  fromType = "ENSEMBL", 
  toType = "ENTREZID", 
  OrgDb = org.Hs.eg.db
)

ekegg <- enrichKEGG(
  gene          = gene_conversion$ENTREZID,
  universe      = universe_conversion$ENTREZID,
  organism      = "hsa", # hsa = Homo sapiens
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

# Convert Entrez ID to Gene Symbol
ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Save KEGG results as table
if(!is.null(ekegg) && nrow(ekegg) > 0) {
  write.table(as.data.frame(ekegg), 
              file = here("results", "differential_expression", "kegg_enrichment_results.tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Visualise KEGG: Dotplot
  kegg_dotplot <- dotplot(ekegg, showCategory = 15) + 
    labs(title = "KEGG Pathway Enrichment")
  ggsave(here("plots", "kegg_dotplot.png"), plot = kegg_dotplot, width = 8, height = 6, dpi = 300)
} else {
  cat("[WARNING] No significant KEGG pathways found.\n")
}


# ==============================================================================
# 6. Gene Set Enrichment Analysis (GSEA) - MSigDB Hallmarks
# ==============================================================================
cat("[INFO] Preparing data for GSEA...\n")

res_df <- res_df %>% 
  mutate(stat_metric = sign(log2FoldChange) * -log10(pvalue)) %>%
  filter(!is.na(stat_metric) & !is.na(clean_id))

# Convert ENSEMBL to ENTREZID for MSigDB
all_conversion <- bitr(res_df$clean_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

gsea_input_df <- res_df %>% 
  inner_join(all_conversion, by = c("clean_id" = "ENSEMBL")) %>%
  arrange(desc(stat_metric))

gene_list <- gsea_input_df$stat_metric
names(gene_list) <- gsea_input_df$ENTREZID

gene_list <- gene_list[!duplicated(names(gene_list))]

cat("[INFO] Fetching MSigDB Hallmark gene sets for Homo sapiens...\n")
# Derive Hallmark (H) categories from из MSigDB database
h_df <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

cat("[INFO] Running GSEA algorithm...\n")
gsea_hallmarks <- GSEA(
  geneList      = gene_list,
  TERM2GENE     = h_df,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  verbose       = FALSE
)

if(!is.null(gsea_hallmarks) && nrow(gsea_hallmarks) > 0) {
  gsea_hallmarks <- setReadable(gsea_hallmarks, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  # Save result as table
  write.table(as.data.frame(gsea_hallmarks), 
              file = here("results", "differential_expression", "gsea_hallmarks_results.tsv"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # --- Visualise GSEA results ---
  
  # 1. Ridgeplot
  gsea_ridge <- ridgeplot(gsea_hallmarks, showCategory = 20) + 
    labs(title = "MSigDB Hallmark Pathways - GSEA Core Enrichment") +
    theme(axis.text.y = element_text(size = 8))
  ggsave(here("plots", "gsea_hallmarks_ridgeplot.png"), plot = gsea_ridge, width = 10, height = 8, dpi = 300)
  
  # 2. GSEA Run Plot (Enrichment Score Profile) for top-5 signal pathways
  top_pathways <- head(gsea_hallmarks@result$ID, 5)
  
  for(i in seq_along(top_pathways)) {
    pathway_name <- top_pathways[i]
    clean_name <- gsub("HALLMARK_", "", pathway_name)
    
    gsea_p <- gseaplot2(gsea_hallmarks, geneSetID = pathway_name, title = pathway_name)
    
    ggsave(here("plots", paste0("gsea_profile_", tolower(clean_name), ".png")), 
           plot = gsea_p, width = 7, height = 6, dpi = 300)
  }
  
  cat("[INFO] GSEA analysis complete! Plots and tables saved successfully.\n")
} else {
  cat("[WARNING] No significant MSigDB Hallmark pathways found via GSEA.\n")
}

cat("[INFO] Functional enrichment analysis completed successfully!\n")
