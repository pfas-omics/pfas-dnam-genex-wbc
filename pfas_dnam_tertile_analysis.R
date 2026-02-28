# ------------------------------------------------------------------------------
# Project: PFAS – DNA Methylation Tertile Analysis
# Date: February 2026
# Script name: pfas_dnam_tertile_analysis.R
# Description: This script performs tertile-based analysis of PFAS exposure 
#              and DNA methylation (450k array), including:
#                - Creating PFAS tertiles
#                - Generating CpG boxplots per PFAS
#                - Linear modeling using limma
#                - Manhattan plots for contrasts
#                - Combining Manhattan plots into panels
#
# Outputs:
#   - Boxplots per PFAS (logit-transformed beta-values)
#   - Manhattan plots for contrasts between PFAS tertiles
#   - 2x2 Manhattan panels for each PFAS
# ------------------------------------------------------------------------------

# ----------------------------
# Install/load packages
# ----------------------------
packages <- c(
  "dplyr", "tibble", "tidyr", "ggplot2", "ggrepel",
  "limma", "patchwork", "png", "grid"
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

prep <- prepare_methylation_matrix("load_your_data.rds")
methylation_matrix <- prep$methylation_matrix
pfas_data <- prep$pfas_data
pfas <- prep$pfas_cols

# ----------------------------
# Create PFAS tertiles
# ----------------------------
pfas_tertiles <- pfas_data %>%
  mutate(across(all_of(pfas), ~ factor(ntile(log(.x), 3), levels = 1:3), .names = "{.col}_tertile"))

pfas_tert <- names(pfas_tertiles)[grep("tertile", names(pfas_tertiles))]

# ----------------------------
# Merge methylation + PFAS tertiles
# ----------------------------
meth_df <- as.data.frame(t(methylation_matrix))
meth_df$labnr <- rownames(meth_df)

final_data_tert <- meth_df %>%
  left_join(pfas_tertiles, by = c("labnr" = "sampleID"))

# ----------------------------
# Boxplot for all CpGs
# ----------------------------

cpg_sign <- data.frame(cpg = colnames(methylation_matrix))  # keep all CpGs
boxplot_df <- final_data_tert %>%
  tidyr::pivot_longer(
    cols = starts_with("cg"),
    names_to = "cpg",
    values_to = "methylation_Value"
  ) %>%
  inner_join(cpg_sign, by = c("cpg"))

# Loop over all PFAS tertiles
for(pfas_to_plot in pfas_tert){
  out_file <- paste0("./PFAS_and_DNAm/boxplot_", pfas_to_plot, ".png")
  png(out_file, width = 2100, height = 1500, res = 600)
  
  ggplot(boxplot_df, aes(x = cpg, y = methylation_Value, fill = .data[[pfas_to_plot]])) +
    geom_boxplot(width = 0.9) +
    labs(y = expression(logit(beta)), x = "", fill = "Tertile") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, color = "black"),
          axis.text.y = element_text(size = 14, color = "black"),
          axis.title.y = element_text(size = 14, color = "black"),
          legend.position = "right",
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15)) +
    scale_fill_brewer(palette = "Dark2") +
    scale_y_continuous(trans = "logit")
  
  dev.off()
  message("Saved boxplot for ", pfas_to_plot)
}
# ----------------------------
# Linear modeling per PFAS tertile
# ----------------------------
pfas_data_tertile <- final_data_tert %>%
  select(c("labnr", pfas_tert, "age", "bmi", "Monocytes","B","CD4T","NK","CD8T","Eosinophils","Neutrophils"))

analyze_pfas_tert <- function(pfas_name, pfas_data, methylation_matrix) {
  rownames(pfas_data) <- pfas_data$labnr
  pfas_data$labnr <- NULL
  
  covariates <- c("age","bmi","Monocytes","B","CD4T","NK","CD8T","Eosinophils","Neutrophils")
  formula_string <- paste0("~ 0 + ", pfas_name, " + ", paste(covariates, collapse = " + "))
  design <- model.matrix(as.formula(formula_string), data = pfas_data)
  
  common_samples <- intersect(colnames(methylation_matrix), rownames(design))
  methylation_matrix_filtered <- methylation_matrix[, common_samples]
  design_filtered <- design[common_samples, ]
  
  fit <- lmFit(methylation_matrix_filtered, design_filtered)
  
  # Define contrasts
  levels <- levels(factor(pfas_data[[pfas_name]]))
  contrast.matrix <- makeContrasts(
    contrasts = c(
      paste0(levels[2], "- ", levels[1]),
      paste0(levels[3], "- ", levels[1]),
      paste0(levels[3], "- ", levels[2])
    ),
    levels = design_filtered
  )
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  top_tables <- lapply(1:ncol(fit2$coefficients), function(i){
    topTable(fit2, coef = i, number = Inf, adjust.method = "fdr", genelist = ann450kSub)
  })
  names(top_tables) <- colnames(fit2$coefficients)
  return(top_tables)
}

results_tert <- lapply(pfas_tert, analyze_pfas_tert, pfas_data = pfas_data_tertile, methylation_matrix = methylation_matrix)
names(results_tert) <- pfas_tert

# ----------------------------
# Manhattan plots
# ----------------------------
plot_manhattan_contrast <- function(pfas_name, contrast_name, tt){
  man <- tt %>%
    select(Name, chr, pos, P.Value, adj.P.Val) %>%
    mutate(logP = -log10(P.Value),
           significant = ifelse(adj.P.Val < 0.05, Name, NA))
  
  man$chr <- gsub("chr","",man$chr)
  man$chr[man$chr=="X"] <- 23
  man$chr[man$chr=="Y"] <- 24
  man$chr <- as.numeric(man$chr)
  
  man <- man[!is.na(man$P.Value) & is.finite(man$P.Value), ]
  man <- man %>% mutate(BP = pos,
                        is_highlight = ifelse(!is.na(significant),"yes","no"),
                        is_annotate = ifelse(adj.P.Val < 0.05,"yes","no"))
  
  don <- man %>%
    group_by(chr) %>% summarise(chr_len=max(BP), .groups="drop") %>%
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    left_join(man, by="chr") %>% arrange(chr,BP) %>%
    mutate(BPcum=BP+tot)
  
  axisdf <- don %>% group_by(chr) %>% summarise(center=(max(BPcum)+min(BPcum))/2)
  
  p <- ggplot(don, aes(BPcum, logP)) +
    geom_point(aes(color=as.factor(chr)), alpha=0.8, size=1.2) +
    scale_color_manual(values=rep(c("grey60","skyblue3"), 12)) +
    scale_x_continuous(breaks=axisdf$center, labels=axisdf$chr) +
    geom_point(data=subset(don,is_highlight=="yes"), color="red", size=2) +
    geom_label_repel(data=subset(don,is_annotate=="yes"), aes(label=Name), size=3.5, max.overlaps=20) +
    coord_cartesian(ylim=c(0,11)) +
    theme_bw() +
    theme(legend.position="none",
          panel.border=element_blank(),
          panel.grid.major.x=element_blank(),
          panel.grid.minor.x=element_blank(),
          axis.title=element_text(size=16),
          axis.text=element_text(size=14)) +
    labs(x="Chromosome", y="-log10(P)", title=paste0(pfas_name," – ",contrast_name))
  
  out_file <- paste0("./PFAS_and_DNAm/Results/manhattan_", pfas_name,"_",contrast_name,".png")
  ggsave(out_file, p, width=14, height=8, dpi=600)
}

# Loop over PFAS tertiles and contrasts
for(pfas_tert in names(results_tert)){
  contrast_names <- names(results_tert[[pfas_tert]])
  for(i in seq_along(contrast_names)){
    plot_manhattan_contrast(pfas_name=pfas_tert,
                            contrast_name=contrast_names[i],
                            tt=results_tert[[pfas_tert]][[i]])
  }
}

# ----------------------------
# Combine Manhattan plots into 2x2 panels
# ----------------------------
png_to_gg <- function(png_file){
  img <- png::readPNG(png_file)
  g <- grid::rasterGrob(img, interpolate=TRUE)
  ggplot() + annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + theme_void()
}

combine_manhattan_panel <- function(pfas_name, file_order, out_dir="./PFAS_and_DNAm/Results/"){
  files <- file.path(out_dir, file_order)
  plots <- lapply(files, png_to_gg)
  combined <- wrap_plots(plots, ncol=2, nrow=2) + plot_annotation(tag_levels="a")
  out_file <- file.path(out_dir, paste0("manhattan_",pfas_name,"_panel.png"))
  ggsave(out_file, combined, width=14, height=12, dpi=600)
  message("Saved 2x2 panel for ", pfas_name)
}
