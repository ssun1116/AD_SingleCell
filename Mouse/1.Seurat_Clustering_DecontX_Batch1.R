rm(list = ls())
library(celda)
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(SingleCellExperiment)

#### 
el01_gex = Read10X(data.dir = "Rawdata/EL01_GEX/cellranger/EL01_GEX_cellranger_count_outs/filtered_feature_bc_matrix/")
el01_sce = SingleCellExperiment(list(counts = el01_gex))
el01_sce = decontX(el01_sce)
el01 <- CreateSeuratObject(round(decontXcounts(el01_sce)), project = "Control_M", min.cells = 10, min.features = 200)
dim(el01) # 20711 9053
el01[["percent.mt"]] <- PercentageFeatureSet(el01, pattern = "^mt-")
el01$log10GenesPerUMI <- log10(el01$nFeature_RNA) / log10(el01$nCount_RNA)
el01 <- subset(el01, subset =  nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5 & log10GenesPerUMI >= 0.8) 
el01 = el01[!row.names(el01) %in% c("Malat1", "Gm26917", "Gm42418", "AY036118")] 

el01 <- NormalizeData(el01, normalization.method = "LogNormalize", scale.factor = 10000)
el01 <- FindVariableFeatures(el01, selection.method = "vst", nfeatures = 2000)
el01 <- ScaleData(el01) # Scale data using all genes
el01 <- RunPCA(el01, npcs = 10, verbose = FALSE)
el01 <- RunUMAP(el01, reduction = "pca", dims = 1:5, min.dist = 0.5)
el01 <- FindNeighbors(el01, reduction = "pca", dims = 1:5)
el01 <- FindClusters(el01, resolution = 0.1)
p1 = DimPlot(el01, reduction="umap", label = T) + NoLegend() 


el02_gex = Read10X(data.dir = "Rawdata/EL02_GEX/cellranger/EL02_GEX_cellranger_count_outs/filtered_feature_bc_matrix/")
el02_sce = SingleCellExperiment(list(counts = el02_gex))
el02_sce = decontX(el02_sce)
el02 <- CreateSeuratObject(round(decontXcounts(el02_sce)), project = "Control_F", min.cells = 10, min.features = 200)
dim(el02)
el02[["percent.mt"]] <- PercentageFeatureSet(el02, pattern = "^mt-")
el02$log10GenesPerUMI <- log10(el02$nFeature_RNA) / log10(el02$nCount_RNA)
el02 <- subset(el02, subset =  nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5 & log10GenesPerUMI >= 0.8) 
el02 = el02[!row.names(el02) %in% c("Malat1", "Gm26917", "Gm42418", "AY036118")] 

el02 <- NormalizeData(el02, normalization.method = "LogNormalize", scale.factor = 10000)
el02 <- FindVariableFeatures(el02, selection.method = "vst", nfeatures = 2000)
el02 <- ScaleData(el02) # Scale data using all genes
el02 <- RunPCA(el02, npcs = 10, verbose = FALSE)
el02 <- RunUMAP(el02, reduction = "pca", dims = 1:5, min.dist = 0.5)
el02 <- FindNeighbors(el02, reduction = "pca", dims = 1:5)
el02 <- FindClusters(el02, resolution = 0.1)
p2 = DimPlot(el02, reduction="umap", label = T) + NoLegend() 


el03_gex = Read10X(data.dir = "Rawdata/EL03_GEX/cellranger/EL03_GEX_cellranger_count_outs/filtered_feature_bc_matrix/")
el03_sce = SingleCellExperiment(list(counts = el03_gex))
el03_sce = decontX(el03_sce)
el03 <- CreateSeuratObject(round(decontXcounts(el03_sce)), project = "Test_M", min.cells = 10, min.features = 200)
dim(el03)
el03[["percent.mt"]] <- PercentageFeatureSet(el03, pattern = "^mt-")
el03$log10GenesPerUMI <- log10(el03$nFeature_RNA) / log10(el03$nCount_RNA)
el03 <- subset(el03, subset =  nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5 & log10GenesPerUMI >= 0.8) 
el03 = el03[!row.names(el03) %in% c("Malat1", "Gm26917", "Gm42418", "AY036118")] 

el03 <- NormalizeData(el03, normalization.method = "LogNormalize", scale.factor = 10000)
el03 <- FindVariableFeatures(el03, selection.method = "vst", nfeatures = 2000)
el03 <- ScaleData(el03) # Scale data using all genes
el03 <- RunPCA(el03, npcs = 10, verbose = FALSE)
el03 <- RunUMAP(el03, reduction = "pca", dims = 1:5, min.dist = 0.5)
el03 <- FindNeighbors(el03, reduction = "pca", dims = 1:5)
el03 <- FindClusters(el03, resolution = 0.1)
p3 = DimPlot(el03, reduction="umap", label = T) + NoLegend() 


el04_gex = Read10X(data.dir = "Rawdata/EL04_GEX/cellranger/EL04_GEX_cellranger_count_outs/filtered_feature_bc_matrix/")
el04_sce = SingleCellExperiment(list(counts = el04_gex))
el04_sce = decontX(el04_sce)
el04 <- CreateSeuratObject(round(decontXcounts(el04_sce)), project = "Test_F", min.cells = 10, min.features = 200)
dim(el04)
el04[["percent.mt"]] <- PercentageFeatureSet(el04, pattern = "^mt-")
el04$log10GenesPerUMI <- log10(el04$nFeature_RNA) / log10(el04$nCount_RNA)
el04 <- subset(el04, subset =  nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5 & log10GenesPerUMI >= 0.8) 
el04 = el04[!row.names(el04) %in% c("Malat1", "Gm26917", "Gm42418", "AY036118")] 

el04 <- NormalizeData(el04, normalization.method = "LogNormalize", scale.factor = 10000)
el04 <- FindVariableFeatures(el04, selection.method = "vst", nfeatures = 2000)
el04 <- ScaleData(el04) # Scale data using all genes
el04 <- RunPCA(el04, npcs = 10, verbose = FALSE)
el04 <- RunUMAP(el04, reduction = "pca", dims = 1:5, min.dist = 0.5)
el04 <- FindNeighbors(el04, reduction = "pca", dims = 1:5)
el04 <- FindClusters(el04, resolution = 0.1)
p4 = DimPlot(el04, reduction="umap", label = T) + NoLegend() 

p = plot_grid(p1, p2, p3, p4, labels = c("Control_Male", "Control_Female", "Test_Male", "Test_Female"))
ggsave(p, file = "Seurat_indv_DecontX_PC.5_min.dist.0.5_res.0.1.pdf", width = 10, height = 10)

seurat.object.list = list()
seurat.object.list[["el01"]] = el01
seurat.object.list[["el02"]] = el02
seurat.object.list[["el03"]] = el03
seurat.object.list[["el04"]] = el04

saveRDS(seurat.object.list, "Seurat_Individual_DecontX.RDS")

