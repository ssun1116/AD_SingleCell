rm(list = ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggrastr)
library(writexl)
library(ggrepel)

combined = readRDS("Data/Seurat_Integrated_DecontX_Edit_DraftAnnot.RDS")
DimPlot(combined, label = T)

## Barplot - Sample proportion
start_color <- "#FFFFFF"  # White
end_color_blue <- "#90C0DF"
end_color_pink <- "#C593C2"

ramp_blue <- colorRampPalette(c(start_color, end_color_blue))
gradient_colors_blue <- ramp_blue(5)
print(gradient_colors_blue)

ramp_pink <- colorRampPalette(c(start_color, end_color_pink))
gradient_colors_pink <- ramp_pink(5)
print(gradient_colors_pink)

p1 = ggplot(combined@meta.data, aes(x=final_annot, fill=cond)) + 
  geom_bar(position = "fill", color = "black") +
  theme_classic(base_size = 10) +
  labs(x = '', y = 'Proportion of cell (%)') +
  scale_fill_manual(values = c("#90C0DF", "#C593C2")) + 
  theme(panel.grid.major = element_blank(),
        axis.text.x=element_text(size = 9, colour="black", hjust=0.5,vjust=0),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.title.y = element_text(size = 10, vjust = 0.2),
        # plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

p2 = ggplot(combined@meta.data, aes(x=final_annot, fill=cond.sex.batch)) + 
  geom_bar(position = "fill", color = "black") +
  theme_classic(base_size = 10) +
  labs(x = '', y = 'Proportion of cell (%)') +
  scale_fill_manual(values = c("#E3EFF7", "#C7DFEF", "#ABCFE7", 
                               "#F0E3EF", "#E2C9E0", "#D3AED1", "#C593C2")) + 
  theme(panel.grid.major = element_blank(),
        axis.text.x=element_text(size = 9, colour="black", hjust=0.5,vjust=0),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.title.y = element_text(size = 10, vjust = 0.2),
        # plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

p = plot_grid(p1, p2, rel_widths = c(0.9, 1))
ggsave(p, file = "Figures/BarPlot_Seurat_Integrated_DraftAnnot_Cluster_Proportion.pdf", width = 15, height = 4)



## DEG Analysis 
combined$sex = sapply(strsplit(as.character(combined$cond.sex),'_'), "[", 2)
combined$batch = sapply(strsplit(as.character(combined$cond.sex.batch),'_'), "[", 3)

set_gene_name <- function(df) {
  df$gene_name = rownames(df)
  return(df)
}


## MAST: Use latent vars
deg.list_mast = list()
for (celltype in unique(combined$final_annot)){
  print(celltype)
  deg = FindMarkers(combined, group.by = "cond", ident.1 = "Test", ident.2 = "Control", 
                    logfc.threshold = -Inf, subset.ident = celltype, 
                    latent.vars = c("sex", "batch"), test.use = "MAST")
  deg.list_mast[[celltype]] = deg
}

# Apply the function to each dataframe in the list
deg.list_mast2 <- lapply(deg.list_mast, set_gene_name)

df_mast = bind_rows(deg.list_mast2, .id = "celltype")
rownames(df_mast) = 1:nrow(df_mast)

saveRDS(df_mast, "Data/DEG_Condition_Mast_Latent_BindList.RDS")


## LR: Use latent vars
deg.list_LR = list()
for (celltype in unique(combined$final_annot)){
  print(celltype)
  deg = FindMarkers(combined, group.by = "cond", ident.1 = "Test", ident.2 = "Control",
                    logfc.threshold = -Inf, subset.ident = celltype,
                    latent.vars = c("sex", "batch"), test.use = "LR")
  deg.list_LR[[celltype]] = deg
}
# Apply the function to each dataframe in the list
deg.list_LR2 <- lapply(deg.list_LR, set_gene_name)

df_LR = bind_rows(deg.list_LR2, .id = "celltype")
rownames(df_LR) = 1:nrow(df_LR)

saveRDS(df_LR, "Data/DEG_Condition_LR_Latent_BindList.RDS")


######################################################
### Visualization

## MAST - latent.vars
for (celltypes in unique(combined$final_annot)){
  
  res1 = df_mast[df_mast$celltype == celltypes,]
  res1$Change = ifelse(res1$p_val_adj < 0.05 & abs(res1$avg_log2FC) > 0.25,
                       ifelse(res1$avg_log2FC > 0, "Up", "Down"), "None")
  x_range = max(abs(res1$avg_log2FC)) * 1.05
  res.gene = res1[res1$Change != "None",]$gene_name
  
  p = ggplot(res1, aes(x = avg_log2FC, y = -log10(p_val_adj) )) + 
    geom_point_rast(aes(fill = Change), 
                    shape = 21, alpha = 0.75,
                    na.rm = F, stroke = 0, size=2.5) + # Make dots bigger
    theme_classic(base_size = 15) + # change theme
    labs(title = paste0('DEG in ', celltypes, ': Test vs. Control (Batch regressed)')) + # Add a title
    xlab(expression(log[2]("Expr A" / "B"))) + # x-axis label
    ylab(expression(-log[10]("FDR"))) + # y-axis label
    geom_hline(yintercept = 1.3, colour = "grey40", linetype='dashed', size = 0.7) + # Add cutoffs
    geom_vline(xintercept = 0.25, linetype='dashed', colour = "grey40", size = 0.7) + # Add 0 lines
    geom_vline(xintercept = -0.25, linetype='dashed', colour = "grey40", size = 0.7) + # Add 0 lines
    ggplot2::xlim(-x_range, x_range) + 
    geom_text_repel(d = subset(res1, res1$gene_name %in% res.gene),
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
  
  ggsave(p, filename = paste0("Figures/DEG_Mast_Latent_Test_vs_Control_", celltypes, ".pdf"),
         width = 6, height = 4.5)
  
}

## LR - latent.vars
for (celltypes in unique(combined$final_annot)){
  
  res1 = df_LR[df_LR$celltype == celltypes,]
  res1$Change = ifelse(res1$p_val_adj < 0.05 & abs(res1$avg_log2FC) > 0.25,
                       ifelse(res1$avg_log2FC > 0, "Up", "Down"), "None")
  x_range = max(abs(res1$avg_log2FC)) * 1.05
  res.gene = res1[res1$Change != "None",]$gene_name
  
  p = ggplot(res1, aes(x = avg_log2FC, y = -log10(p_val_adj) )) + 
    geom_point_rast(aes(fill = Change), 
                    shape = 21, alpha = 0.75,
                    na.rm = F, stroke = 0, size=2.5) + # Make dots bigger
    theme_classic(base_size = 15) + # change theme
    labs(title = paste0('DEG in ', celltypes, ': Test vs. Control (Batch regressed)')) + # Add a title
    xlab(expression(log[2]("Expr A" / "B"))) + # x-axis label
    ylab(expression(-log[10]("FDR"))) + # y-axis label
    geom_hline(yintercept = 1.3, colour = "grey40", linetype='dashed', size = 0.7) + # Add cutoffs
    geom_vline(xintercept = 0.25, linetype='dashed', colour = "grey40", size = 0.7) + # Add 0 lines
    geom_vline(xintercept = -0.25, linetype='dashed', colour = "grey40", size = 0.7) + # Add 0 lines
    ggplot2::xlim(-x_range, x_range) + 
    geom_text_repel(d = subset(res1, res1$gene_name %in% res.gene),
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
  
  ggsave(p, filename = paste0("Figures/DEG_LR_Latent_Test_vs_Control_", celltypes, ".pdf"),
         width = 6, height = 4.5)
  
}

