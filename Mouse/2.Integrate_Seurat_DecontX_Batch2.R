rm(list = ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
library(cowplot)

fs.list = readRDS('Data/SeuratList_Individual_DecontX_batch2_0209.RDS')

# ## Basic merge
# merged <- Reduce(merge, fs.list)
# merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
# merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(merged) # Scale data using all genes
# merged <- ScaleData(merged, features = all.genes)
# merged <- RunPCA(merged, features = VariableFeatures(object = merged))
# merged <- merged %>% RunUMAP(reduction = "pca", dims = 1:5, min.dist = 0.3) %>% 
#   FindNeighbors(dims = 1:5) %>% 
#   FindClusters(resolution = 0.4)
# 
# ## Harmony
# merged2 = Reduce(merge, fs.list)
# merged2 <- NormalizeData(merged2, normalization.method = "LogNormalize", scale.factor = 10000)
# merged2 <- FindVariableFeatures(merged2, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(merged2) # Scale data using all genes
# merged2 <- ScaleData(merged2, features = all.genes)
# merged2 <- RunPCA(merged2, features = VariableFeatures(object = merged2))
# merged2 <- merged2 %>% RunHarmony("orig.ident", plot_convergence = TRUE) %>% 
#   RunUMAP(reduction = "harmony", dims = 1:5, min.dist = 0.3) %>% 
#   FindNeighbors(reduction = "harmony", dims = 1:5) %>% 
#   FindClusters(resolution = 0.4) 
# 
# DimPlot(object = merged2, reduction = "pca", pt.size = .1, group.by = "orig.ident")

## Integrate 

# ## Sampling cells - for identical number of cells per sample
# fs.list2 = list()
# sample.size = 6000
# set.seed(111)
# 
# for (file in names(fs.list)){
#   obj = fs.list[[file]]
#   orig.size = ncol(obj) 
#   sampled.cells <- colnames(obj)[sample(x = orig.size, size = sample.size, replace = F)]
#   obj.sub <- subset(obj, cells = sampled.cells)
#   fs.list2[[file]] = obj.sub
# }

fs.list2 = fs.list

features <- SelectIntegrationFeatures(object.list = fs.list2)
anchors <- FindIntegrationAnchors(object.list = fs.list2, anchor.features = features)
combined <- IntegrateData(anchorset = anchors)
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 10, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:5)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:5)
combined <- FindClusters(combined, resolution = 0.5)

saveRDS(combined, "Data/Seurat_Combined_DecontX_batch2_0209.RDS")

## DimPlot

p1 = DimPlot(combined, label = TRUE, pt.size = .1)
p2 = DimPlot(combined, split.by = "orig.ident", label = TRUE, pt.size = .1)

ggsave(p1, file = "Figures/DimPlot_Seurat_Combined_DecontX_batch2.pdf", width = 6, height = 5)
ggsave(p2, file = "Figures/DimPlot_Seurat_Combined_DecontX_batch2_Split.pdf", width = 12.5, height = 5)



