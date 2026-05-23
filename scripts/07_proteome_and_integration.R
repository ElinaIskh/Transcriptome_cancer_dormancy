# ==============================================================================
# Project: Integrative Transcriptomic and Proteomic Analysis of Cancer Dormancy
# Script: 07_proteome_and_integration.R
# Description: Multi-week proteomic analysis (cytoplasm & nucleus) and dynamic
#              integration with RNA-Seq (DESeq2) data.
# ==============================================================================

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(readxl)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)
library(here)

# Create output folders if they do not exist
mkdir_if_missing <- function(path) if(!dir.exists(path)) dir.create(path, recursive = TRUE)
mkdir_if_missing(here("results", "proteom"))
mkdir_if_missing(here("plots"))

# Normalisation function
normalize_counts <- function(x) { (x / sum(x, na.rm = TRUE)) * 1e6 }

# ==============================================================================
# 1. Dowload and preprocessing of mass-spectrometry data
# ==============================================================================
cat("[INFO] Loading and preprocessing proteomic data...\n")

# Download data
prot_excel_path <- here("results", "proteom", "proteome_results.xlsx")
if(!file.exists(prot_excel_path)) stop("Error: 'proteome_results.xlsx' missing in 'proteom/'!")

cyt_raw <- read_excel(prot_excel_path, sheet = "Cytoplasm")
nucl_raw <- read_excel(prot_excel_path, sheet = "Nucleus")

cyt_col_name <- colnames(cyt_raw)[1]
nucl_col_name <- colnames(nucl_raw)[1]

# Cytoplasm fraction
dea_results_cyt <- cyt_raw %>%
  mutate(Gene = str_extract(!!sym(cyt_col_name), "(?<=GN=)[^ ]+")) %>%
  filter(!is.na(Gene)) %>%
  dplyr::select(Gene, Accession, starts_with("cyt")) %>%
  mutate(across(where(is.numeric), normalize_counts)) %>%
  mutate(
    lfc_1w = log2((cyt_1w + 1) / (cyt_ctrl + 1)),
    lfc_2w = log2((cyt_2w + 1) / (cyt_ctrl + 1)),
    lfc_3w = log2((cyt_3w + 1) / (cyt_ctrl + 1))
  )

# Nucleus fraction
dea_results_nucl <- nucl_raw %>%
  mutate(Gene = str_extract(!!sym(nucl_col_name), "(?<=GN=)[^ ]+")) %>%
  filter(!is.na(Gene)) %>%
  dplyr::select(Gene, Accession, starts_with("nucl")) %>%
  mutate(across(where(is.numeric), normalize_counts)) %>%
  mutate(
    lfc_1w = log2((nucl_1w + 1) / (nucl_ctrl + 1)),
    lfc_2w = log2((nucl_2w + 1) / (nucl_ctrl + 1)),
    lfc_3w = log2((nucl_3w + 1) / (nucl_ctrl + 1))
  )

# ==============================================================================
# 2. Dowload transcriptome data and integration
# ==============================================================================
cat("[INFO] Integrating with RNA-Seq results...\n")

rna_path <- here("results", "differential_expression", "differential_expression_results.tsv")
if(!file.exists(rna_path)) stop("Error: Run script 04 (DESeq2) first!")
rna_df <- read.table(rna_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Функция для интеграции конкретной фракции протеома с РНК
integrate_omics <- function(prot_df, fraction_name) {
  inner_join(
    rna_df %>% dplyr::select(Gene = gene_name, log2FC_RNA = log2FoldChange, padj_RNA = padj),
    prot_df %>% dplyr::select(Gene, Accession, lfc_1w, lfc_2w, lfc_3w),
    by = "Gene"
  ) %>% mutate(Fraction = fraction_name)
}

integrated_cyt <- integrate_omics(dea_results_cyt, "Cytoplasm")
integrated_nucl <- integrate_omics(dea_results_nucl, "Nucleus")
all_integrated <- bind_rows(integrated_cyt, integrated_nucl)

write.table(all_integrated, here("results", "proteom", "transcriptome_proteome_timecourse.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)


# ==============================================================================
# 4. Visualisation: Multi-quadrant graph
# ==============================================================================
cat("[INFO] Generating cross-omics scatter plots...\n")

plt_df <- integrated_cyt %>%
  mutate(Group = case_when(
    log2FC_RNA > 1 & lfc_3w > 1 ~ "Concordant Up",
    log2FC_RNA < -1 & lfc_3w < -1 ~ "Concordant Down",
    abs(log2FC_RNA) > 1 | abs(lfc_3w) > 1 ~ "Discordant / Delayed",
    TRUE ~ "Stable"
  ))

top_labels <- plt_df %>% filter(Group != "Stable") %>% arrange(desc(abs(log2FC_RNA) + abs(lfc_3w))) %>% head(10)

scatter_p <- ggplot(plt_df, aes(x = log2FC_RNA, y = lfc_3w, color = Group)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_text_repel(data = top_labels, aes(label = Gene), color = "black", size = 3) +
  scale_color_manual(values = c("Concordant Up"="#d73027", "Concordant Down"="#4575b4", "Discordant / Delayed"="#fe9929", "Stable"="#grey80")) +
  labs(title = "RNA-Seq vs Cytoplasm Proteome (Week 3)", x = "RNA log2FC (mTOR treated)", y = "Protein log2FC (Week 3 / Control)") +
  theme_bw()

ggsave(here("plots", "rna_vs_protein_week3_cyt.png"), plot = scatter_p, width = 7, height = 6, dpi = 300)

# ==============================================================================
# 5. ComplexHeatmap
# ==============================================================================
cat("[INFO] Building multi-week ComplexHeatmap...\n")

heatmap_subset <- integrated_cyt %>%
  filter(padj_RNA < 0.05) %>%
  arrange(desc(abs(log2FC_RNA))) %>%
  head(40)

if(nrow(heatmap_subset) > 0) {
  mat_rna <- as.matrix(heatmap_subset$log2FC_RNA)
  mat_prot <- as.matrix(heatmap_subset[, c("lfc_1w", "lfc_2w", "lfc_3w")])
  rownames(mat_rna) <- rownames(mat_prot) <- heatmap_subset$Gene
  colnames(mat_prot) <- c("Prot 1w", "Prot 2w", "Prot 3w")
  
  col_fun <- colorRamp2(c(-3, 0, 3), c("#4575b4", "white", "#d73027"))
  
  ht_rna <- Heatmap(mat_rna, name = "RNA log2FC", col = col_fun, cluster_columns = FALSE, column_title = "Transcriptome")
  ht_prot <- Heatmap(mat_prot, name = "Protein log2FC", col = col_fun, cluster_columns = FALSE, column_title = "Proteome Dynamic")
  
  png(here("plots", "integrated_timecourse_heatmap.png"), width = 2000, height = 2400, res = 300)
  draw(ht_rna + ht_prot, main_title = "Dynamic Shift of Transcriptome vs Cytoplasm Proteome")
  dev.off()
}
