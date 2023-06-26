rm(list = ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
library(cowplot)

fs.list = readRDS('Seurat_Individual_SoupX_0622.RDS')

## Sampling cells - for identical number of cells per sample
fs.list2 = list()
sample.size = 7000
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

DimPlot(object = merged, reduction = "umap", pt.size = .1, group.by = "orig.ident")
saveRDS(merged, "Data/data_harmony_SoupX_Sampled_0622.RDS")


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

DimPlot(object = combined, reduction = "umap", pt.size = .1, group.by = "orig.ident")

#saveRDS(combined, "Data/data_combined_SoupX_Sampled_0622.RDS")


## DimPlot
p1 = DimPlot(combined, group.by = "orig.ident", pt.size = .1)
p2 = DimPlot(merged, group.by = "orig.ident", pt.size = .1)

p = plot_grid(p1, p2, ncol = 2)

ggsave(p, file = "Seurat_integrate_SoupX.pdf", width = 12, height = 5)


