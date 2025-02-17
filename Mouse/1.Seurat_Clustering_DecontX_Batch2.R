rm(list = ls())
library(celda)
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(SingleCellExperiment)
setwd("~/Dropbox/ANGPT2_AD_v2/")

#### 
EL05_gex = Read10X(data.dir = "Data/Rawdata/EL05/cellranger/EL05_cellranger_count_outs/filtered_feature_bc_matrix/")
EL05_sce = SingleCellExperiment(list(counts = EL05_gex))
EL05_sce = decontX(EL05_sce)
EL05 <- CreateSeuratObject(round(decontXcounts(EL05_sce)), project = "Control_M", min.cells = 10, min.features = 200)
dim(EL05) # 20711 9053
EL05[["percent.mt"]] <- PercentageFeatureSet(EL05, pattern = "^mt-")
EL05$log10GenesPerUMI <- log10(EL05$nFeature_RNA) / log10(EL05$nCount_RNA)
VlnPlot(EL05, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

EL05 <- subset(EL05, subset =  nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5 & log10GenesPerUMI >= 0.8) 
EL05 = EL05[!row.names(EL05) %in% c("Malat1", "Gm26917", "Gm42418", "AY036118")] 

EL05 <- NormalizeData(EL05, normalization.method = "LogNormalize", scale.factor = 10000)
EL05 <- FindVariableFeatures(EL05, selection.method = "vst", nfeatures = 2000)
EL05 <- ScaleData(EL05) # Scale data using all genes
EL05 <- RunPCA(EL05, npcs = 10, verbose = FALSE)
EL05 <- RunUMAP(EL05, reduction = "pca", dims = 1:5, min.dist = 0.5)
EL05 <- FindNeighbors(EL05, reduction = "pca", dims = 1:5)
EL05 <- FindClusters(EL05, resolution = 0.1)
p1 = DimPlot(EL05, reduction="umap", label = T) + NoLegend() 


EL07_gex = Read10X(data.dir = "Data/Rawdata/EL07/cellranger/EL07_cellranger_count_outs/filtered_feature_bc_matrix/")
EL07_sce = SingleCellExperiment(list(counts = EL07_gex))
EL07_sce = decontX(EL07_sce)
EL07 <- CreateSeuratObject(round(decontXcounts(EL07_sce)), project = "Test_M", min.cells = 10, min.features = 200)
dim(EL07)
EL07[["percent.mt"]] <- PercentageFeatureSet(EL07, pattern = "^mt-")
EL07$log10GenesPerUMI <- log10(EL07$nFeature_RNA) / log10(EL07$nCount_RNA)
VlnPlot(EL07, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))


EL07 <- subset(EL07, subset =  nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5 & log10GenesPerUMI >= 0.8) 
EL07 = EL07[!row.names(EL07) %in% c("Malat1", "Gm26917", "Gm42418", "AY036118")] 

EL07 <- NormalizeData(EL07, normalization.method = "LogNormalize", scale.factor = 10000)
EL07 <- FindVariableFeatures(EL07, selection.method = "vst", nfeatures = 2000)
EL07 <- ScaleData(EL07) # Scale data using all genes
EL07 <- RunPCA(EL07, npcs = 10, verbose = FALSE)
EL07 <- RunUMAP(EL07, reduction = "pca", dims = 1:5, min.dist = 0.5)
EL07 <- FindNeighbors(EL07, reduction = "pca", dims = 1:5)
EL07 <- FindClusters(EL07, resolution = 0.1)
p2 = DimPlot(EL07, reduction="umap", label = T) + NoLegend() 


EL08_gex = Read10X(data.dir = "Data/Rawdata/EL08/cellranger/EL08_cellranger_count_outs/filtered_feature_bc_matrix/")
EL08_sce = SingleCellExperiment(list(counts = EL08_gex))
EL08_sce = decontX(EL08_sce)
EL08 <- CreateSeuratObject(round(decontXcounts(EL08_sce)), project = "Test_F", min.cells = 10, min.features = 200)
dim(EL08)
EL08[["percent.mt"]] <- PercentageFeatureSet(EL08, pattern = "^mt-")
EL08$log10GenesPerUMI <- log10(EL08$nFeature_RNA) / log10(EL08$nCount_RNA)
VlnPlot(EL08, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))


EL08 <- subset(EL08, subset =  nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5 & log10GenesPerUMI >= 0.8) 
EL08 = EL08[!row.names(EL08) %in% c("Malat1", "Gm26917", "Gm42418", "AY036118")] 

EL08 <- NormalizeData(EL08, normalization.method = "LogNormalize", scale.factor = 10000)
EL08 <- FindVariableFeatures(EL08, selection.method = "vst", nfeatures = 2000)
EL08 <- ScaleData(EL08) # Scale data using all genes
EL08 <- RunPCA(EL08, npcs = 10, verbose = FALSE)
EL08 <- RunUMAP(EL08, reduction = "pca", dims = 1:5, min.dist = 0.5)
EL08 <- FindNeighbors(EL08, reduction = "pca", dims = 1:5)
EL08 <- FindClusters(EL08, resolution = 0.1)
p3 = DimPlot(EL08, reduction="umap", label = T) + NoLegend() 

p = plot_grid(p1, p2, p3, labels = c("Control_Male", "Test_Male", "Test_Female"))
ggsave(p, file = "Figures/DimPlot_Seurat_Individual_DecontX_PC.5_mindist.0.5_res.0.1_batch2.pdf", width = 10, height = 10)

seurat.object.list = list()
seurat.object.list[["EL05"]] = EL05
seurat.object.list[["EL07"]] = EL07
seurat.object.list[["EL08"]] = EL08

saveRDS(seurat.object.list, "Data/SeuratList_Individual_DecontX_batch2_0209.RDS")

