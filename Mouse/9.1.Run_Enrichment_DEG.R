rm(list = ls())
library(readxl)
library(dplyr)
library(annotate)
library(ggplot2)
source("~run_enrichment_v8_GOBP.R")

df_mast = readRDS("Data/DEG_Condition_Mast_Latent_BindList.RDS")
colnames(df_mast) 
# [1] "celltype"   "p_val"      "avg_log2FC" "pct.1"      "pct.2"      "p_val_adj"  "gene_name" 

## cluster
for (cluster in unique(df_mast$celltype)) {
  tryCatch({
    print(cluster)
    deg_cluster_tmp = df_mast[df_mast$celltype == cluster,]
    res = run_all(deg_cluster_tmp,
                  organism = 'mouse',
                  lfc_column = "avg_log2FC",
                  padj_column = "p_val_adj",
                  gene_name_column = 'gene_name',
                  lfc_thresh = 0.25,
                  padj_cutoff = 0.05,
                  minGSSize = 100,
                  maxGSSize = 1000,
                  output_file = paste0('Tables/cluster_', cluster, '_enrichment_result.xlsx'))

    res_pl = plot_all(res,
                      width = 20,
                      height = 12.5,
                      lfc_thresh = 0.25,
                      padj_cutoff = 0.05,
                      output_file = paste0('Figures/cluster_', cluster, '_enrichment_result.pdf'))
  }, error = function(e) {
    cat("An error occurred for cluster", cluster, ":", conditionMessage(e), "\n")
  })
}
