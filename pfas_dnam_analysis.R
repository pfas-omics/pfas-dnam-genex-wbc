# ------------------------------------------------------------------------------
# Project: PFAS – DNA Methylation Analysis Workflow
# Date: February 2026
# Script name: pfas_dnam_analysis.R
# Description: This script performs DNA methylation analysis in relation to PFAS 
#              exposure using Illumina 450k data. Key steps include:
#                - Preparing methylation matrix and PFAS data
#                - Linear modeling of PFAS vs CpG methylation (limma)
#                - Manhattan and volcano plot generation
#                - CpG overlap analysis across PFAS
#                - Gene Ontology (GO) enrichment analysis (missMethyl)
#                - UpSetR plots for shared CpGs across PFAS
#
# Outputs:
#   - Methylation vs PFAS linear model results
#   - Manhattan plots per PFAS
#   - Volcano plots per PFAS
#   - UpSet plot for overlapping significant CpGs
#   - GO enrichment tables for significant CpGs
# ------------------------------------------------------------------------------

# ----------------------------
# Install/load packages
# ----------------------------

packages <- c(
  "dplyr", "tibble", "limma", "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylation450kmanifest", "ggplot2", "ggrepel", 
  "tidyr", "UpSetR", "missMethyl"
)

for(pkg in packages){
  if(!requireNamespace(pkg, quietly = TRUE)){
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# ----------------------------
# Load annotation
# ----------------------------
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k %>%
  data.frame() %>%
  select("chr", "pos", "strand", "Name", "Probe_rs", "Probe_maf", "CpG_rs", "CpG_maf", "SBE_rs", "SBE_maf",                 
         "Islands_Name", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Accession",  
         "UCSC_RefGene_Group", "Phantom", "DMR", "Enhancer", "HMM_Island", "Regulatory_Feature_Name",
          "Regulatory_Feature_Group", "DHS")

anno_df <- ann450kSub %>%
  as.data.frame() %>%
  rownames_to_column("Name") %>%
  select(Name, chr, pos)

# ----------------------------
# Prepare methylation matrix + PFAS data
# ----------------------------
prepare_methylation_matrix <- function(
  dnam_data_path,
  sample_id_col = "sampleID",
  covariates = c("age","bmi","smoking_status","Monocytes","B","CD4T","NK","CD8T","Eosinophils","Neutrophils"),
  pfas_pattern = "^PF",
  select_cols = NULL
) {
  dnam_data <- readRDS(dnam_data_path)
  
  if(!is.null(select_cols)){
    dnam_subset <- dnam_data %>% select(all_of(c(sample_id_col, select_cols)))
  } else {
    dnam_subset <- dnam_data
  }
  
  dnam_subset <- dnam_subset %>% filter(complete.cases(across(all_of(covariates))))
  
  pfas_cols <- names(dnam_subset)[grep(pfas_pattern, names(dnam_subset))]
  pfas_data <- dnam_subset %>% select(all_of(c(sample_id_col, covariates, pfas_cols)))
  final_labnr <- dnam_subset[[sample_id_col]]
  
  methylation_matrix <- dnam_data %>%
    filter(.data[[sample_id_col]] %in% final_labnr) %>%
    select(all_of(sample_id_col), starts_with("cg")) %>%
    t() %>% as.data.frame()
  
  colnames(methylation_matrix) <- methylation_matrix[1, ]
  methylation_matrix <- methylation_matrix[-1, ]
  methylation_matrix <- methylation_matrix[, colnames(methylation_matrix) %in% final_labnr]
  colnames(methylation_matrix) <- final_labnr
  methylation_matrix <- as.matrix(methylation_matrix)
  
  return(list(
    methylation_matrix = methylation_matrix,
    pfas_data = pfas_data,
    pfas_cols = pfas_cols,
    final_labnr = final_labnr
  ))
}

# ----------------------------
# PFAS linear model analysis
# ----------------------------
analyze_pfas_linear <- function(
  pfas_name,
  pfas_data,
  methylation_matrix,
  sample_id_col = "sampleID",
  covariates = NULL,
  top_n = 20,
  pval_cutoff = 0.001
) {
  aligned_idx <- match(colnames(methylation_matrix), pfas_data[[sample_id_col]])
  aligned_pfas_data <- pfas_data[aligned_idx, ]
  
  if(!is.null(covariates)){
    keep <- complete.cases(aligned_pfas_data[, c(pfas_name, covariates), drop=FALSE])
  } else {
    keep <- complete.cases(aligned_pfas_data[, pfas_name, drop=FALSE])
  }
  
  aligned_pfas_data <- aligned_pfas_data[keep, ]
  methylation_matrix <- methylation_matrix[, keep]
  
  formula_string <- if(!is.null(covariates)){
    paste0("~", pfas_name, " + ", paste(covariates, collapse=" + "))
  } else {
    paste0("~", pfas_name)
  }
  
  design <- model.matrix(as.formula(formula_string), data = aligned_pfas_data)
  
  fit <- lmFit(methylation_matrix, design)
  fit <- eBayes(fit)
  
  all_dmps <- topTable(fit, coef = pfas_name, adjust.method = "BH", number = Inf, sort.by = "P")
  top_dmps <- topTable(fit, coef = pfas_name, adjust.method = "BH", number = top_n, sort.by = "P")
  filtered_dmps <- all_dmps[all_dmps$P.Value < pval_cutoff, ]
  
  return(list(All_DMP = all_dmps, Top_DMPs = top_dmps, Filtered_DMPs = filtered_dmps))
}

# ----------------------------
# Prepare and run analysis
# ----------------------------
prep <- prepare_methylation_matrix("load_your_data.rds")
methylation_matrix <- prep$methylation_matrix
pfas_data <- prep$pfas_data
pfas <- prep$pfas_cols

covariates_no_smoking <- c("age","bmi","Monocytes","B","CD4T","NK","CD8T","Neutrophils")
covariates_with_smoking <- c(covariates_no_smoking,"smoking_status")

results_no_smoking <- lapply(pfas,
  analyze_pfas_linear,
  pfas_data = pfas_data,
  methylation_matrix = methylation_matrix,
  sample_id_col = "sampleID",
  covariates = covariates_no_smoking
)
names(results_no_smoking) <- pfas

# Save outputs
# ===============================
# saveRDS(results_no_smoking, "pfas_DMPS_results.rds")
# ------------------------------------------------------------------------------

results_with_smoking <- lapply(pfas,
  analyze_pfas_linear,
  pfas_data = pfas_data,
  methylation_matrix = methylation_matrix,
  sample_id_col = "sampleID",
  covariates = covariates_with_smoking
)
names(results_with_smoking) <- pfas


# ----------------------------
# Annotate and plot Manhattan
# ----------------------------
results_no_smoking_annotated <- lapply(results_no_smoking, function(x){
  x$All_DMP <- x$All_DMP %>%
    tibble::rownames_to_column("Name") %>%
    left_join(anno_df, by="Name")
  x
})

out_dir <- "./Output/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

plot_manhattan_element <- function(element_name){
  pfas_manhattan <- results_no_smoking_annotated[[element_name]]$All_DMP
  pfas_man <- pfas_manhattan %>%
    select(Name, chr, pos, P.Value, adj.P.Val) %>%
    mutate(logP = -log10(P.Value),
           significant = ifelse(adj.P.Val < 0.05, Name, NA))
  
  # Chromosome numeric
  pfas_man$chr <- gsub("chr","",pfas_man$chr)
  pfas_man$chr[pfas_man$chr=="X"] <- 23
  pfas_man$chr[pfas_man$chr=="Y"] <- 24
  pfas_man$chr <- as.numeric(pfas_man$chr)
  pfas_man <- pfas_man[!is.na(pfas_man$P.Value) & is.finite(pfas_man$P.Value), ]
  
  pfas_man <- pfas_man %>%
    mutate(BP = pos,
           is_highlight = ifelse(!is.na(significant),"yes","no"),
           is_annotate = ifelse(adj.P.Val < 0.05,"yes","no"))
  
  don <- pfas_man %>%
    group_by(chr) %>%
    summarise(chr_len = max(BP)) %>%
    mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
    select(-chr_len) %>%
    left_join(pfas_man, ., by=c("chr"="chr")) %>%
    arrange(chr,BP) %>%
    mutate(BPcum = BP + tot)
  
  axisdf <- don %>% group_by(chr) %>% summarise(center=(max(BPcum)+min(BPcum))/2)
  n_chr <- length(unique(don$chr))
  
  manhattan_plot <- ggplot(don, aes(x=BPcum, y=logP)) +
    geom_point(aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
    scale_color_manual(values=rep(c("grey","skyblue"), length.out=n_chr)) +
    scale_x_continuous(label=axisdf$chr, breaks=axisdf$center) +
    scale_y_continuous(expand=c(0,0)) +
    geom_point(data=subset(don,is_highlight=="yes"), color="red", size=2) +
    {if(nrow(subset(don,is_annotate=="yes"))>0) geom_label_repel(
      data=subset(don,is_annotate=="yes"),
      aes(label=Name), size=4, max.overlaps=20)} +
    theme_bw() +
    theme(legend.position="none",
          panel.border=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          axis.title=element_text(size=16),
          axis.text=element_text(size=14)) +
    labs(x="Chromosome", y="-log10(P)") +
    coord_cartesian(ylim=c(0,10.5))
  
  out_file <- paste0(out_dir,"manhattan_",element_name,".png")
  ggsave(out_file, manhattan_plot, width=14, height=8, dpi=600)
  message("Saved Manhattan plot for ", element_name)
}

elements <- names(results_no_smoking_annotated)
lapply(elements, plot_manhattan_element)

# ----------------------------
# Volcano plots
# ----------------------------
plot_volcano <- function(dmp_results, output_file, logFC_limits=c(-0.25,0.25), pval_cutoff=0.05){
  volcano_data <- dmp_results$All_DMP %>%
    tibble::rownames_to_column("cpg") %>%
    select(cpg, logFC, P.Value, adj.P.Val, UCSC_RefGene_Name) %>%
    mutate(significant = ifelse(adj.P.Val < pval_cutoff, "Significant", "Not Significant"),
           negLogP = -log10(P.Value)) %>%
    filter(!is.na(significant))
  
  png(output_file, width=2100, height=1500, res=600)
  
  p <- ggplot(volcano_data, aes(x=logFC, y=negLogP, color=significant)) +
    geom_point(size=3) +
    labs(x=expression(Log[2] ~ "Fold Change"),
         y=expression(-log[10] ~ "("~italic(p)~")"),
         color="") +
    xlim(logFC_limits) +
    scale_color_manual(values=c("Significant"="red","Not Significant"="grey")) +
    theme_minimal() +
    theme(axis.text=element_text(size=14,color="black"),
          axis.title=element_text(size=16,color="black"),
          legend.title=element_text(size=14),
          legend.text=element_text(size=12),
          legend.position="right") +
    geom_hline(yintercept=-log10(pval_cutoff), linetype="dashed", color="black") +
    geom_vline(xintercept=c(-0.01,0.01), linetype="dashed", color="black")
  
  print(p)
  dev.off()
  message("Saved volcano plot to ", output_file)
}

for(pfas_name in names(results_no_smoking)){
  output_file <- paste0(out_dir,"volcano_",pfas_name,".png")
  plot_volcano(results_no_smoking[[pfas_name]], output_file)
}

# ----------------------------
# CpG overlap & GO enrichment
# ----------------------------
pfas_names <- names(results_dmps)
pval_cutoff <- 0.01

filtered_cpgs_list <- lapply(pfas_names, function(pfas){
  df <- results_dmps[[pfas]]$All_DMP
  df_filtered <- df %>% filter(P.Value < pval_cutoff) %>%
    mutate(PFAS = pfas) %>%
    tibble::rownames_to_column("Name") %>%
    select(Name, PFAS)
  df_filtered
})
names(filtered_cpgs_list) <- pfas_names

combined_pval <- bind_rows(filtered_cpgs_list) %>%
  distinct() %>%
  filter(complete.cases(.))

combined_pval_binary <- combined_pval %>%
  mutate(value=1) %>%
  pivot_wider(names_from=PFAS, values_from=value, values_fill=list(value=0))
rownames(combined_pval_binary) <- combined_pval_binary$Name
combined_pval_binary <- combined_pval_binary %>% select(-Name)

pdf(file=paste0(out_dir,"upset_filtered_cpgs.pdf"), width=18, height=14)
upset(as.data.frame(combined_pval_binary),
      order.by="degree",
      nsets=length(pfas_names),
      line.size=0,
      text.scale=4.5,
      point.size=8,
      sets.bar.color="black",
      matrix.color="black",
      main.bar.color="black",
      sets=rev(pfas_names),
      keep.order=TRUE)
dev.off()

overlap <- combined_pval %>% group_by(Name) %>%
  summarise(compound_count = n_distinct(PFAS)) %>% arrange(desc(compound_count))

# GO enrichment
go_results_list <- lapply(pfas_names, function(pfas){
  sig_cpgs <- filtered_cpgs_list[[pfas]]$Name
  gometh(sig.cpg=sig_cpgs, collection="GO") %>%
    tibble::rownames_to_column("ID") %>%
    mutate(PFAS=pfas)
})
names(go_results_list) <- pfas_names

all_go <- do.call(rbind, go_results_list)

# GO for all overlapping CpGs
all_overlap_cpgs <- unique(combined_pval$Name)
go_all_overlap <- gometh(sig.cpg=all_overlap_cpgs, sig.genes=TRUE, collection="GO") %>%
  tibble::rownames_to_column("ID")

top_n <- 20
top_go_by_pfas <- all_go %>% filter(P.DE < 0.01) %>%
  group_by(PFAS) %>%
  arrange(P.DE) %>%
  slice_head(n=top_n)

write.table(top_go_by_pfas, file=paste0(out_dir,"go_terms_top20.txt"),
            quote=FALSE, row.names=FALSE, sep="\t")
