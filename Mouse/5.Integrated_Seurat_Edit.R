rm(list = ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
library(wesanderson)
library(cowplot)

combined = readRDS("Data/Seurat_Integrated_DecontX_AllBatch.RDS")
DimPlot(combined, label = T)

## Add 0.7 resolution
combined <- FindClusters(combined, resolution = 0.7)

# Apply is.na() to each column in combined@meta.data
result <- apply(combined@meta.data, MARGIN = 2, FUN = function(x) sum(is.na(x)))
print(result)

combined@meta.data = combined@meta.data[,c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "log10GenesPerUMI", "seurat_clusters", "barcode", "celltype", "orig.ident_batch", "batch")]

combined@meta.data$cond = ifelse(combined@meta.data$orig.ident %in% c("Control_M", "Control_F"), "Control", "Test")
combined$cond.sex = combined$orig.ident
combined$cond.sex = factor(combined$cond.sex, levels = c("Control_M", "Control_F", "Test_M", "Test_F"))
combined$cond.sex.batch = combined$orig.ident_batch
combined$cond.sex.batch = factor(combined$cond.sex.batch, levels = c("Control_M_Batch1","Control_M_Batch2", "Control_F_Batch1", "Test_M_Batch1",  "Test_M_Batch2", "Test_F_Batch1", "Test_F_Batch2"))

table(combined$cond) # Condition
table(combined$cond.sex) # Condition+Sex
table(combined$cond.sex.batch) # Condition+Sex+Batch

combined$orig.ident = ifelse(combined$cond.sex.batch == "Control_M_Batch1", "EL01",
                             ifelse(combined$cond.sex.batch == "Control_M_Batch2", "EL05", 
                                    ifelse(combined$cond.sex.batch == "Control_F_Batch1", "EL02",
                                           ifelse(combined$cond.sex.batch == "Test_M_Batch1", "EL03", 
                                                  ifelse(combined$cond.sex.batch == "Test_M_Batch2", "EL07",
                                                         ifelse(combined$cond.sex.batch == "Test_F_Batch1", "EL04", "EL08"))))))

table(combined$orig.ident)

combined$orig.ident_batch <- NULL
combined$batch <- NULL
combined$barcode <- NULL

## Barplot visualization
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
        axis.text.x=element_text(size = 9, colour="black", angle = 45, hjust=1,vjust=1),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.title.y = element_text(size = 10, vjust = 0.2),
        # plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

ggsave(p1, file = "Figures/ANGPT2_AD/BarPlot_Seurat_Integrated_Cluster_Proportion_Condition_0827.pdf",
       width = 6, height = 4)


p2 = ggplot(combined@meta.data, aes(x=seurat_clusters, fill=cond.sex.batch)) + 
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
ggsave(p, file = "Figures/BarPlot_Seurat_Integrated_Cluster_Proportion.pdf",
       width = 12, height = 4)

combined$celltype = ifelse(combined$celltype == "DG", "Ex.DG", combined$celltype)

p1 = DimPlot(combined, group.by = "seurat_clusters", label = TRUE)
p2 = DimPlot(combined, group.by = "celltype", label = TRUE)

p = plot_grid(p1, p2)
ggsave(p, file = "Figures/DimPlot_Seurat_Integrated_DecontX_AllBatch_SingleR.pdf", width = 12.5, height = 5)

p_split = DimPlot(combined, split.by = "celltype", ncol = 5, label = TRUE)
ggsave(p_split, file = "Figures/DimPlot_Seurat_Integrated_DecontX_AllBatch_SingleR_Split.pdf",
       width = 12.5, height = 7.5)

saveRDS(combined, "Data/Seurat_Integrated_DecontX_AllBatch_Edit.RDS")


