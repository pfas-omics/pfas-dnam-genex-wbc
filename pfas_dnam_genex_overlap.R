# -------------------------------------------------------------------------------
# Project: PFAS – DNAm & Gene Expression Overlap
# Date: February 2026
# Script: pfas_dnam_genex_overlap.R
# Description: Generalized workflow to analyze overlap between DNAm and gene expression associated with PFAS exposures. 
# Outputs: 
#   - Table of overlapping CpG–gene–PFAS pairs
#   - Spearman correlation plot of all overlapping CpG–gene pairs per PFAS.
# -------------------------------------------------------------------------------

required_packages <- c("tidyverse", "ggplot2")

for(pkg in required_packages){
  if(!requireNamespace(pkg, quietly = TRUE)){
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}
# ----------------------------
# Load previously saved results
# ----------------------------
dmps_results <- readRDS("dmps_results.rds")
degs_results <- readRDS("degs_results.rds")

# ----------------------------
# Load methylation and expression matrices
# ----------------------------
# Make sure you have saved these during the previous steps
methylation_matrix <- readRDS("methylation_matrix.rds")  # CpGs x samples
expression_matrix <- readRDS("expression_matrix.rds")    # genes x samples

# ----------------------------
# Flatten DMPs to CpG-gene pairs
# ----------------------------
all_dmps <- bind_rows(dmps_results, .id="PFAS") %>%
  filter(!is.na(UCSC_RefGene_Name)) %>%
  separate_rows(UCSC_RefGene_Name, sep=";") %>%
  mutate(gene = trimws(UCSC_RefGene_Name)) %>%
  select(PFAS, Name, gene)

# Flatten DEGs
all_degs <- bind_rows(degs_results, .id="PFAS") %>%
  select(PFAS, gene)

# ----------------------------
# Identify overlapping PFAS-gene-CpG pairs
# ----------------------------
overlap_df <- all_dmps %>%
  inner_join(all_degs, by=c("PFAS","gene"))

# ----------------------------
# Subset matrices to overlapping CpGs and genes
# ----------------------------
cpgs_to_keep <- unique(overlap_df$Name)
genes_to_keep <- unique(overlap_df$gene)

meth_sub <- methylation_matrix[rownames(methylation_matrix) %in% cpgs_to_keep, , drop=FALSE]
expr_sub <- expression_matrix[rownames(expression_matrix) %in% genes_to_keep, , drop=FALSE]

# ----------------------------
# Compute Spearman correlations for each pair
# ----------------------------
df_pairs <- overlap_df %>%
  filter(Name %in% rownames(meth_sub) & gene %in% rownames(expr_sub)) %>%
  distinct(Name, gene, PFAS)

results_correlation_pairs <- df_pairs %>%
  rowwise() %>%
  mutate(
    rho = cor(as.numeric(meth_sub[Name, ]), as.numeric(expr_sub[gene, ]), use="complete.obs", method="spearman"),
    pval = cor.test(as.numeric(meth_sub[Name, ]), as.numeric(expr_sub[gene, ]), use="complete.obs", method="spearman")$p.value
  ) %>%
  ungroup() %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  mutate(sig_label = ifelse(padj < 0.05, "p < 0.05", "p ≥ 0.05"),
         pair = paste(Name, gene, sep="-"),
         pair_italic = paste0("italic('", pair, "')"))

# ----------------------------
# Plot CpG–gene correlations
# ----------------------------
p <- ggplot(results_correlation_pairs, aes(x = reorder(pair, rho), y = rho, color = sig_label)) +
  geom_point(size = 6) +
  coord_flip() +
  labs(x = "", y = "Spearman correlation", color = "") +
  scale_color_manual(values = c("p ≥ 0.05" = "cornflowerblue", "p < 0.05" = "brown2")) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x  = element_text(size = 14, colour = "black", margin = margin(t = 25)),
    axis.text.y  = element_text(size = 14, colour = "black")
  ) +
  scale_x_discrete(labels = parse(text = results_correlation_pairs$pair_italic))

# Save plot
ggsave(
  filename = "CpG_gene_correlations.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 600
)

# ----------------------------
# Save overlapping table
# ----------------------------
write.table(overlap_df, "CpG_DEG_overlap_table.txt", sep="\t", row.names=FALSE)

message("Overlap analysis complete. Correlation plot saved as 'CpG_gene_correlations.png'.")
