## Pathway testing
rm(list = ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(readxl)
source("run_enrichment_v8_GOBP.R") # GOBP

############## ROSMAP - Endothelial ###############

filename = "DESeq2_Pseudobulk_Endothelial_Diagnosis_cov.age_death.pmi"

df = read_excel(paste0("Tables/ROSMAP/", filename, ".xlsx"))

res = run_all(df, 
              organism = 'human', 
              lfc_column = "log2FoldChange", 
              padj_column = "pvalue", 
              gene_name_column = 'gene_name', 
              lfc_thresh = 0.25, 
              padj_cutoff = 0.05, 
              minGSSize = 100, 
              maxGSSize = 1000,
              output_file = paste0('Tables/ROSMAP/Enrichment_', filename, '_240812.xlsx'))

res_pl = plot_all(res,
                  width = 15,
                  height = 10,
                  lfc_thresh = 0.25, 
                  padj_cutoff = 0.05, 
                  output_file = paste0('Figures/ROSMAP/Enrichment_', filename, '_240521.pdf'))

res_data = res$res_list$gsea_go_bp

p1 = gseaplot(res_data, geneSetID = 14, by = "runningScore", 
         title = res_data$Description[14])

p2 = gseaplot2(res_data, geneSetID = 90,
         title = res_data$Description[90])

p3 = gseaplot(res_data, geneSetID = 34, by = "runningScore", 
         title = res_data$Description[34])

p4 = gseaplot(res_data, geneSetID = 39, by = "runningScore", 
         title = res_data$Description[39])

p5 = gseaplot(res_data, geneSetID = 47, by = "runningScore", 
         title = res_data$Description[47])

p6 = gseaplot(res_data, geneSetID = 67, by = "runningScore", 
         title = res_data$Description[67])

p7 = gseaplot(res_data, geneSetID = 85, by = "runningScore", 
         title = res_data$Description[85])

p8 = gseaplot(res_data, geneSetID = 94, by = "runningScore", 
         title = res_data$Description[94])

p9 = gseaplot(res_data, geneSetID = 122, by = "runningScore", 
         title = res_data$Description[122])

p10 = gseaplot(res_data, geneSetID = 123, by = "runningScore", 
         title = res_data$Description[123])

p11 = gseaplot(res_data, geneSetID = 129, by = "runningScore", 
         title = res_data$Description[129])

plot_grid(p1, p2, p3, p4, p5, p6, p7, NA, p8, p9, p10, p11, ncol = 4)
ggsave(last_plot(), filename = "GSEA_EnrichmentPlot_Endo_DEG_240503.pdf",
       width = 4, height = 2)



############## ROSMAP - DV4 ANG2 High vs. Low ###############

filename = "DESeq2_Pseudobulk_Endothelial_Diagnosis_cov.age_death.pmi"

fs <- list.files('Tables/ROSMAP/DESeq2_Pseudo_ANG2_Group_AD/', recursive = T)
fs1 = fs[grepl("ANG2_High_vs_Low", fs)]
fs1 = gsub(".xlsx", "", fs1)

for (filename in fs1){
  
  df = read_excel(paste0("Tables/ROSMAP/DESeq2_Pseudo_ANG2_Group_AD/", filename, ".xlsx"))
  
  tryCatch({
    res = run_all(df, 
                  organism = 'human', 
                  lfc_column = "log2FoldChange", 
                  padj_column = "pvalue", 
                  gene_name_column = 'gene_name', 
                  lfc_thresh = 0.25, 
                  padj_cutoff = 0.05, 
                  minGSSize = 100, 
                  maxGSSize = 1000,
                  output_file = paste0('Tables/ROSMAP/Enrichment_', filename, '.xlsx'))
  
    res_pl = plot_all(res,
                      width = 15,
                      height = 10,
                      lfc_thresh = 0.25, 
                      padj_cutoff = 0.05, 
                      output_file = paste0('Figures/ROSMAP/Enrichment_', filename, '.pdf'))
    }, error = function(e) {
      cat("An error occurred for cluster", filename, ":", conditionMessage(e), "\n")
    })
}


# ## cluster
# for (cluster in unique(deg_cluster$cluster)) {
#   tryCatch({
#     deg_cluster_tmp = deg_cluster[deg_cluster$cluster == cluster,]
#     res = run_all(deg_cluster_tmp, 
#                   organism = 'mouse', 
#                   lfc_column = "avg_log2FC", 
#                   padj_column = "p_val_adj", 
#                   gene_name_column = 'gene', 
#                   lfc_thresh = 0.2, 
#                   padj_cutoff = 0.05, 
#                   minGSSize = 100, 
#                   maxGSSize = 1000,
#                   output_file = paste0('cluster_', cluster, '_enrichment_result.xlsx'))
#     
#     res_pl = plot_all(res,
#                       width = 20,
#                       height = 17.5,
#                       lfc_thresh = 0.2, 
#                       padj_cutoff = 0.05, 
#                       output_file = paste0('cluster_', cluster, '_enrichment_result.pdf'))
#   }, error = function(e) {
#     cat("An error occurred for cluster", cluster, ":", conditionMessage(e), "\n")
#   })
# }
# 
# 
# ## deg_list
# for (name in names(deg_list)){
#   tryCatch({
#     deg_list_tmp = deg_list[[name]]
#     deg_list_tmp$gene_name = rownames(deg_list_tmp)
#     res = run_all(deg_list_tmp, 
#                   organism = 'mouse', 
#                   lfc_column = "avg_log2FC", 
#                   padj_column = "p_val_adj", 
#                   gene_name_column = 'gene_name', 
#                   lfc_thresh = 0.2, 
#                   padj_cutoff = 0.05, 
#                   minGSSize = 100, 
#                   maxGSSize = 1000,
#                   output_file = paste0('list_', name, '_enrichment_result.xlsx'))
#     
#     res_pl = plot_all(res,
#                       width = 20,
#                       height = 17.5,
#                       lfc_thresh = 0.2, 
#                       padj_cutoff = 0.05, 
#                       output_file = paste0('list_', name, '_enrichment_result.pdf'))
#   }, error = function(e) {
#     cat("An error occurred for cluster", name, ":", conditionMessage(e), "\n")
#   })
# }
# 
