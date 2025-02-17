rm(list = ls())
library(readxl)
library(dplyr)
library(annotate)
library(ggplot2)
setwd("~/Dropbox/ANGPT2_AD_v2")
source("~/Dropbox/CodeShareHub/EnrichmentPlot/run_enrichment_v4.3_GOBP.R")

df_mast = readRDS("Data/DEG_Condition_Mast_Latent_BindList_240219.RDS")
colnames(df_mast) 
# [1] "celltype"   "p_val"      "avg_log2FC" "pct.1"      "pct.2"      "p_val_adj"  "gene_name" 

for (celltypes in unique(df_mast$celltype)){
  print("===== Run analysis =====")
  print(celltypes)
  res = df_mast[df_mast$celltype == celltypes,]

  ## 1. Plot volcano
  print("1. Plot volcano")
  res$Change = ifelse(res$p_val_adj < 0.05 & abs(res$avg_log2FC) > 0.25,
                      ifelse(res$avg_log2FC > 0, "Up", "Down"), "None")
  x_range = max(abs(res$avg_log2FC)) * 1.05
  res.gene = res[res$Change != "None",]$gene_name
  
  p = ggplot(res, aes(x = avg_log2FC, y = -log10(p_val_adj) )) + 
    geom_point_rast(aes(fill = Change), 
                    shape = 21, alpha = 0.75,
                    na.rm = F, stroke = 0, size=2.5) + # Make dots bigger
    theme_classic(base_size = 15) + # change theme
    labs(title = paste0('DEG in ', celltypes, ': Test vs. Control (Batch Regressed)')) + # Add a title
    xlab(expression(log[2]("Expr A" / "B"))) + # x-axis label
    ylab(expression(-log[10]("FDR"))) + # y-axis label
    geom_hline(yintercept = 1.3, colour = "grey40", linetype='dashed', size = 0.7) + # Add cutoffs
    geom_vline(xintercept = 0.25, linetype='dashed', colour = "grey40", size = 0.7) + # Add 0 lines
    geom_vline(xintercept = -0.25, linetype='dashed', colour = "grey40", size = 0.7) + # Add 0 lines
    ggplot2::xlim(-x_range, x_range) + 
    geom_text_repel(d = subset(res, res$gene_name %in% res.gene),
                    aes(label = gene_name),
                    size = 3.5) +
    theme(
      panel.grid.major = element_line(colour = "grey80", size = 0.5),
      panel.grid.minor = element_line(colour = "grey90", size = 0.25),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.title = element_text(size =12),
      title = element_text(size = 12),
      #legend.position = c(0.5, 0.9),
      #legend.direction = "horizontal",
      legend.background = element_rect(fill = "white"),
      legend.margin=margin(0.1,0.1,0.1,0.1)) +
    scale_fill_manual(values = c("Up" = "#E64B35", 
                                 "Down" = "#3182bd", 
                                 "None" = "grey"))
  
  print(paste0("Number of DEGs - Up: ", nrow(res[res$Change == "Up",]),
               " Down: ", nrow(res[res$Change == "Down",])))
  
  
  ## 2. GSEA
  print("2. Run GSEA")
  
  ## feature 1: numeric vector
  geneList = res[,"avg_log2FC"]
  ## feature 2: named vector
  names(geneList) = as.character(res$gene_name)
  ## feature 3: decreasing order
  geneList = sort(geneList, decreasing = TRUE)
  
  ego_gse <- gseGO(geneList     = geneList,
                   OrgDb        = org.Mm.eg.db,
                   keyType = "SYMBOL",
                   ont          = "BP",
                   minGSSize    = 100,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   eps = 0,
                   verbose      = FALSE)
  ego_gse_df = ego_gse@result %>% filter(p.adjust < 0.05)
  print(paste0("Number of Enriched GSEA terms: ", nrow(ego_gse_df)))
  
  # ## 3. ORA 
  # print("3. Run ORA")
  # 
  # ## ORA input
  # res.up = res[res$p_val_adj < 0.05 & res$avg_log2FC > 0.25,]$gene_name
  # res.dn = res[res$p_val_adj < 0.05 & res$avg_log2FC < -0.25,]$gene_name
  # 
  # ## ORA - Up
  # ego_up <- enrichGO(gene          = res.up,
  #                    OrgDb         = org.Mm.eg.db,
  #                    keyType = "SYMBOL",
  #                    ont           = "BP",
  #                    pAdjustMethod = "bonferroni",
  #                    pvalueCutoff  = 0.01,
  #                    qvalueCutoff  = 0.05,
  #                    readable      = TRUE)
  # ego_up_df = ego_up@result %>% filter(p.adjust < 0.05)
  # print(paste0("Number of Enriched ORA-Up terms: ", nrow(ego_up_df)))
  # 
  # ## ORA - Dn
  # ego_dn <- enrichGO(gene          = res.dn,
  #                    OrgDb         = org.Mm.eg.db,
  #                    keyType = "SYMBOL",
  #                    ont           = "BP",
  #                    pAdjustMethod = "bonferroni",
  #                    pvalueCutoff  = 0.01,
  #                    qvalueCutoff  = 0.05,
  #                    readable      = TRUE)
  # ego_dn_df = ego_dn@result %>% filter(p.adjust < 0.05)
  # print(paste0("Number of Enriched ORA-Down terms: ", nrow(ego_dn_df)))
  
  ## 4. Plot enrichemt results
  print("4. Plot enrichment results") 
  df_GSEA = ego_gse@result
  # df_GSEA$Description_ID = paste(df_GSEA$Description, " (", df_GSEA$ID, ")", sep = "")
  df_GSEA = df_GSEA[df_GSEA$p.adjust < 0.05,]
  df_GSEA = df_GSEA %>% top_n(20, -log10(p.adjust))
  df_GSEA$binary_Dir = ifelse(df_GSEA$NES > 0, 1, 0)
  df_GSEA$Description = ifelse(df_GSEA$NES > 0, paste0(df_GSEA$Description, "   "),
                                      paste0("   ", df_GSEA$Description))

  p_gsea = ggplot(df_GSEA, aes(x = NES, y = reorder(Description, NES))) +
    geom_bar(stat = "identity", 
             color = 'black', lwd = 0.5,
             aes(fill = NES), alpha= 1,
             position = "identity", show.legend = F) +
    scale_fill_gradient2(low = '#169194', mid = 'white', high = '#C593C2', midpoint = 0) + 
    geom_text(aes(y=Description, x=0, label= Description), hjust=df_GSEA$binary_Dir) + 
    xlim(c(-max(abs(df_GSEA$NES)), max(abs(df_GSEA$NES)))) +
    #coord_fixed(ratio = 0.15) +
    labs(x = 'NES', y = '') +
    #theme_minimal(base_size = 7) +
    geom_vline(xintercept = 0, colour = "black", size = 0.5) + # Add 0 lines
    theme_classic(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      axis.text.x = element_text(size = 9, colour="black", vjust = 1, hjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 12, vjust = 0.2),
      axis.line.y = element_blank(),
      plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
      plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm"),
      #plot.margin = margin(t = 10, b = 10)
      )
  
  # tryCatch({
  #   edox_up <- pairwise_termsim(ego_up)
  #   if (dim(edox_up)[1] > 5){
  #     p1 <- treeplot(edox_up)
  #   } else {
  #     p1 <- ggplot() + theme_void()
  #   }
  # },
  # error = function(e) {
  #   # Code to handle the error
  #   p1 <- ggplot() + theme_void()
  # })
  # 
  # tryCatch({
  #   edox_dn <- pairwise_termsim(ego_dn)
  #   if (dim(edox_dn)[1] > 5){
  #     p2 <- treeplot(edox_dn)
  #   } else {
  #     p2 <- ggplot() + theme_void()
  #   }
  # },
  # error = function(e) {
  #   # Code to handle the error
  #   p2 <- ggplot() + theme_void()
  # })
  
  
  ## Merge plots
  p_upper = plot_grid(p, p_gsea, labels = c("", "GSEA"), nrow = 1, rel_widths = c(1, 1.5))
  #p_lower = plot_grid(p1, p2, labels = c("ORA-Up", "ORA-Down"), ncol = 1)
  #p_final = plot_grid(p_upper, p_lower, ncol = 1, rel_heights = c(1, 1.75))
  
  ## Save
  ggsave(p_upper, filename = paste0("Figures/Enrichment_Figure_", celltypes,"_multipanel.pdf"),
         width = 14, height = 4)
  
  ego_list = list()
  # ego_list[["ORA_Up"]] = ego_up
  # ego_list[["ORA_Down"]] = ego_dn
  ego_list[["GSEA"]] = ego_gse
  
  res_list = list(DEGs = res, 
                  # ORA_Up = ego_up@result, 
                  # ORA_Down = ego_dn@result, 
                  GSEA = ego_gse@result)
  
  write_xlsx(x = res_list, path = paste0("Tables/Enrichment_Table_", celltypes, "_multipanel.xlsx"))
  
  print("===== Complete =====")
  print("")
}


# 
# df = res[["GSEA_ego"]]@result
# View(df)
# rownames(df) = 1:252
# 
# p1 = gseaplot(res[["GSEA_ego"]], by = "runningScore", geneSetID = 178, title = firstup(res[["GSEA_ego"]]$Description[178]))
# 
# p2 = gseaplot(res[["GSEA_ego"]], by = "runningScore", geneSetID = 207, title = firstup(res[["GSEA_ego"]]$Description[207]))
# 
# plot_grid(p1, p2)
# 
#
