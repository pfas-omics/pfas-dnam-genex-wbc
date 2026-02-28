# PFAS – DNA Methylation and Gene Expression Analysis


## Project Overview
This project investigates the relationship between **PFAS exposure, DNA methylation, and gene expression** using linear models and enrichment analyses.

It includes:
- Analysis of **DNA CpG methylation** changes associated with **PFAS**
- **Gene expression** profiling associated with **PFAS**
- Overlap analysis between **DNA methylation** and **gene expression**
- Examine associations between PFAS exposure and **estimated immune cell proportions**, both to investigate alterations in circulating immune cells and to contextualize observed methylation changes


## Project Structure

1. **DNA methylation (DNAm) analysis**
   - Script: **`pfas_dnam_analysis.R`**
   - Linear models of PFAS vs CpG methylation using **limma**.
   - Generates **Manhattan** and **Volcano plots**.
   - Conducts **GO enrichment analysis** for significant CpGs.
   - Outputs:
     - `pfas_DMPS_results.rds` → All DMP results per PFAS.
     - `manhattan_*.png` → Manhattan plots per PFAS.
     - `volcano_*.png` → Volcano plots per PFAS.
     - `go_terms_top20.txt` → Top GO terms per PFAS.

2. **DNAm analysis Tertile Analysis**
   - Script: **`pfas_dnam_tertile_analysis.R`**
   - DNA methylation CpGs are analyzed by **PFAS exposure tertiles**.
   - Linear models using **limma** compare methylation between **low, medium, and high PFAS tertiles**.
   - Generates **Volcano plots** and **summary tables** per tertile comparison.
   - Outputs:
     - `pfas_DMPS_tertiles_results.rds` → DMPs per PFAS tertile.
     - `volcano_tertiles_*.png` → Volcano plots for tertile comparisons.
  
3. **Gene expression Analysis**
   - Script: **pfas_gene_expression.R**
   - Linear models using **limma** with covariates including **age, BMI, and WBC proportions**.
   - Outputs:
     - `pfas_expression_results.rds` → All DEGs per PFAS.
     - `filtered_DEGs_all_PFAS.txt` → Filtered DEGs (p < 0.05) across all PFAS.

4. **Overlapp between DNAm and gene expression analyses**
   - Script: **`pfas_dnam_genex_overlap.R`**
   - Matches **DNAm and Gene expression results** to identify overlapping CpG–gene pairs.
   - Computes **Spearman correlations** between methylation and expression for overlapping pairs.
   - Generates plots of correlations.
   - Outputs:
     - `CpG_DEG_overlap_table.txt` → Table of overlapping CpG–gene–PFAS pairs.
     - `CpG_gene_correlations.png` → Spearman correlation plot of overlapping pairs.
    
5. **Associations between PFAS and WBC**
   - Linear models of **PFAS concentrations vs estimated immune cell proportions**.
   - Adjusted for **age, BMI**, and other covariates.
   - Outputs:
     - `pfas_wbc_results.rds` → Regression coefficients and p-values for PFAS vs WBC.
  

## Usage

1. **Input data**:
   - RDS files, containing:
        - DNAm + PFAS + covariates.
        - gene expression + PFAS + covariates.

