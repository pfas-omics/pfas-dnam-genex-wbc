# ------------------------------------------------------------------------------
# Project: PFAS and Gene Expression Analysis
# Date: February 2026
# Script name: pfas_gene_expression.R
# Description:
#   Generalized workflow for PFAS–gene expression association analysis using limma, including:
#      - Run expression models
#      - Combine DEGs for all PFAS
# Outputs:
#   - Top DEGs per PFAS
#   - Filtered DEGs by p-value
# ------------------------------------------------------------------------------

# ===============================
# Load / install packages
# ===============================
required_packages <- c("tidyverse", "limma")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# ===============================
# Define parameters
# ===============================

# Input expression RDS file
expression_rds <- "./GitHub/expression_Data.rds"

# Sample ID column
sample_id_col <- "sampleID"

# Columns containing expression values
expression_range <- c("AADACL1.1", "ZZZ3")  # first to last gene column

# PFAS column pattern (regex)
pfas_pattern <- "^PF"

# Covariates for linear models
covariates <- c(
  "age", "bmi",
  "Monocytes", "B", "CD4T", "NK", "CD8T", "Neutrophils"
)

# DEG filtering
pval_cutoff <- 0.05
top_n <- 20

# ===============================
# Load & prepare expression data
# ===============================
final_exp <- readRDS(expression_rds)

# Expression matrix (genes × samples)
expression_matrix <- final_exp %>%
  select(
    all_of(sample_id_col),
    all_of(expression_range[1]):all_of(expression_range[2])
  ) %>%
  column_to_rownames(sample_id_col) %>%
  t() %>%
  as.matrix()

# Phenotype data (samples × covariates + PFAS)
pheno_data <- final_exp %>%
  select(
    all_of(sample_id_col),
    matches(pfas_pattern),
    all_of(covariates)
  )

# List of PFAS
pfas <- names(pheno_data)[grepl(pfas_pattern, names(pheno_data))]

# ===============================
# General PFAS → expression model function
# ===============================
run_expression_model <- function(
  pfas_name,
  pheno_data,
  expression_matrix,
  sample_id_col,
  covariates,
  pval_cutoff = 0.05,
  top_n = 20
) {
  # Align samples
  idx <- match(colnames(expression_matrix), pheno_data[[sample_id_col]])
  aligned_pheno <- pheno_data[idx, ]
  
  # Drop incomplete cases
  keep <- complete.cases(aligned_pheno[, c(pfas_name, covariates)])
  aligned_pheno <- aligned_pheno[keep, ]
  exp_mat <- expression_matrix[, keep]
  
  # Design matrix
  design <- model.matrix(
    as.formula(
      paste("~", pfas_name, "+", paste(covariates, collapse = " + "))
    ),
    data = aligned_pheno
  )
  
  # Fit linear model with empirical Bayes
  fit <- lmFit(exp_mat, design)
  fit <- eBayes(fit)
  
  # Extract results
  all_degs <- topTable(
    fit,
    coef = pfas_name,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  ) %>%
    rownames_to_column("gene")
  
  list(
    All_DEGs = all_degs,
    Top_DEGs = slice_head(all_degs, n = top_n),
    Filtered_DEGs = filter(all_degs, P.Value < pval_cutoff)
  )
}

# ===============================
# Run analysis for all PFAS
# ===============================
expression_results <- lapply(
  pfas,
  run_expression_model,
  pheno_data = pheno_data,
  expression_matrix = expression_matrix,
  sample_id_col = sample_id_col,
  covariates = covariates,
  pval_cutoff = pval_cutoff,
  top_n = top_n
)

names(expression_results) <- pfas

# ===============================
# Combine filtered DEGs for all PFAS
# ===============================
filtered_degs_all <- map_dfr(
  names(expression_results),
  function(pf) {
    expression_results[[pf]]$Filtered_DEGs %>%
      mutate(PFAS = pf)
  }
)

# ===============================
# Save outputs
# ===============================
# saveRDS(expression_results, "pfas_expression_results.rds")
# write.table(filtered_degs_all, "filtered_DEGs_all_PFAS.txt", sep = "\t", row.names = FALSE)
# ------------------------------------------------------------------------------
