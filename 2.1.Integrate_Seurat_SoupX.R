rm(list = ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
library(cowplot)

fs.list = readRDS('Data/Seurat_Individual_SoupX_Doublet_0626.RDS')

## Sampling cells - for identical number of cells per sample
fs.list2 = list()
sample.size = 7500
set.seed(111)
for (file in c("el01", "el02", "el03", "el04")){
  obj = fs.list[[file]]
  orig.size = ncol(obj) 
  sampled.cells <- colnames(obj)[sample(x = orig.size, size = sample.size, replace = F)]
  obj.sub <- subset(obj, cells = sampled.cells)
  fs.list2[[file]] = obj.sub
}

# Harmony
merged = Reduce(merge, fs.list2)
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(merged) # Scale data using all genes
merged <- ScaleData(merged, features = all.genes)
merged <- RunPCA(merged, features = VariableFeatures(object = merged))
merged <- merged %>% RunHarmony("orig.ident", plot_convergence = TRUE) %>%
  RunUMAP(reduction = "harmony", dims = 1:5, min.dist = 0.3) %>%
  FindNeighbors(reduction = "harmony", dims = 1:5) %>%
  FindClusters(resolution = 0.4)

merged@meta.data$Doublet = coalesce(merged@meta.data$DF.classifications_0.25_21_601, 
                                    merged@meta.data$DF.classifications_0.25_2_510, 
                                    merged@meta.data$DF.classifications_0.25_2_498, 
                                    merged@meta.data$DF.classifications_0.25_32_612)
table(merged@meta.data$Doublet)
p1 = DimPlot(object = merged, reduction = "umap", pt.size = .1, group.by = "Doublet")
saveRDS(merged, "Data/data_harmony_SoupX_Doublet_Sampled_0626.RDS")


## Integrate 
features <- SelectIntegrationFeatures(object.list = fs.list2)
anchors <- FindIntegrationAnchors(object.list = fs.list2, anchor.features = features)
combined <- IntegrateData(anchorset = anchors)
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- FindNeighbors(combined, dims = 1:5)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, dims = 1:5)

combined@meta.data$Doublet = coalesce(combined@meta.data$DF.classifications_0.25_21_601, 
                                      combined@meta.data$DF.classifications_0.25_2_510, 
                                      combined@meta.data$DF.classifications_0.25_2_498, 
                                      combined@meta.data$DF.classifications_0.25_32_612)
table(combined@meta.data$Doublet)
p2 = DimPlot(object = combined, reduction = "umap", pt.size = .1, group.by = "Doublet")

DefaultAssay(combined) = "RNA"
saveRDS(combined, "Data/data_combined_SoupX_Doublet_Sampled_0626.RDS")


## DimPlot
p1 = DimPlot(combined, group.by = "orig.ident", pt.size = .1)
p2 = DimPlot(merged, group.by = "orig.ident", pt.size = .1)
p = plot_grid(p1, p2, ncol = 2)

ggsave(p, file = "Seurat_integrate_SoupX.pdf", width = 12, height = 5)


