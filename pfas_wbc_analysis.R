# ------------------------------------------------------------------------------
# Project: PFAS – WBC Analysis
# Date: February 2026
# Script name: pfas_wbc_analysis.R
# Description: Generalized workflow to analyze associations between immune cell proportions 
# and PFAS exposures. 
# Outputs: linear model results for each PFAS-cell type combination.
# ------------------------------------------------------------------------------

packages <- c(
  "tidyverse", "dplyr", "tibble", "readr", "ggplot2", "tidyr", "purrr"
)

for(pkg in packages){
  if(!requireNamespace(pkg, quietly = TRUE)){
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# ----------------------------
# 1. Define PFAS and immune cell columns
# ----------------------------
pfas_vars <- names(pfas_data)[grep("^PF", names(pfas_data))]        # PFAS exposures
immune_cells <- c("Monocytes","B","CD4T","NK","CD8T","Eosinophils","Neutrophils")    # change if needed

# ----------------------------
# Generate all combinations of PFAS and immune cell types
# ----------------------------
combinations <- expand.grid(pfas = pfas_vars,
                            cell_type = immune_cells,
                            stringsAsFactors = FALSE)

# ----------------------------
# Function to run linear model
# ----------------------------
run_model <- function(pfas, cell_type, covariates = c("age", "bmi"), data) {
  formula_string <- as.formula(paste(cell_type, "~", pfas, "+", paste(covariates, collapse = " + ")))
  model <- lm(formula_string, data = data)
  coefs <- summary(model)$coefficients
  
  # Return tidy results
  df <- tibble(
    pfas = pfas,
    cell_type = cell_type,
    term = rownames(coefs),
    estimate = coefs[, "Estimate"],
    std_error = coefs[, "Std. Error"],
    t_value = coefs[, "t value"],
    p_value = coefs[, "Pr(>|t|)"]
  )
  return(df)
}

# ----------------------------
# Run models for all combinations
# ----------------------------
results_list <- purrr::pmap(combinations, run_model, data = pfas_data)

# Combine all results into a single data frame
all_results <- bind_rows(results_list)

# ----------------------------
# Save or view results
# ----------------------------
write_csv(all_results, "Output/pfas_vs_immune_cells.csv")
print(all_results)
