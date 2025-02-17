rm(list = ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
library(cowplot)

combined_batch1 = readRDS("Data/Seurat_Combined_DecontX_Annotated_final_batch1_0212.RDS")
combined_batch1$batch = "Batch1"
combined_batch1$orig.ident_batch = paste(combined_batch1$orig.ident, combined_batch1$batch, sep = "_")

combined_batch2 = readRDS('Data/Seurat_Combined_DecontX_Annotated_final_batch2_0212.RDS')
combined_batch2$batch = "Batch2"
combined_batch2$orig.ident_batch = paste(combined_batch2$orig.ident, combined_batch2$batch, sep = "_")

combined_all = merge(combined_batch1, combined_batch2) 
table(combined_all$orig.ident)
table(combined_all$orig.ident_batch)

# split the dataset into a list of two seurat objects (stim and CTRL)
combined.list <- SplitObject(combined_all, split.by = "orig.ident")

# Integrate
features <- SelectIntegrationFeatures(object.list = combined.list)
anchors <- FindIntegrationAnchors(object.list = combined.list, anchor.features = features)
combined <- IntegrateData(anchorset = anchors)
DefaultAssay(combined) <- "integrated"

# npcs = 
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 5, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:5)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:5)
combined <- FindClusters(combined, resolution = 0.7)
DimPlot(combined, label = TRUE, pt.size = .1)

# combined <- RunPCA(combined, npcs = 10, verbose = FALSE)
# combined <- RunUMAP(combined, reduction = "pca", dims = 1:5)
# combined <- FindNeighbors(combined, reduction = "pca", dims = 1:5)
# combined <- FindClusters(combined, resolution = 0.5)

saveRDS(combined, "Data/Seurat_Integrated_DecontX_AllBatch_0212.RDS")

## DimPlot
p1 = DimPlot(combined, label = TRUE, pt.size = .1)

combined@meta.data$orig.ident_batch = factor(combined@meta.data$orig.ident_batch, levels = c("Control_M_Batch1", "Control_F_Batch1", "Test_M_Batch1", "Test_F_Batch1",
                                                                                                "Control_M_Batch2", "Test_M_Batch2", "Test_F_Batch2"))
p2 = DimPlot(combined, split.by = "orig.ident_batch", ncol = 4, 
             label = TRUE, pt.size = .1,)

ggsave(p1, file = "Figures/DimPlot_Seurat_Integrated_DecontX_AllBatch.pdf", width = 6.5, height = 5)
ggsave(p2, file = "Figures/DimPlot_Seurat_Integrated_DecontX_AllBatch_Split.pdf", width = 12.5, height = 7)



