rm(list = ls())
library(celda)
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(SoupX)
library(SingleCellExperiment)
setwd("~/Dropbox/ANGPT2_AD/")

##### EL01 #####
## Prepare data for SoupX
el01_toc = Read10X(data.dir = "Rawdata/EL01_GEX/cellranger/EL01_GEX_cellranger_count_outs/filtered_feature_bc_matrix/")
el01_tod = Read10X(data.dir = "Rawdata/EL01_GEX/cellranger/EL01_GEX_cellranger_count_outs/raw_feature_bc_matrix/")
el01_sc = SoupChannel(el01_tod, el01_toc)
el01_sc

## Prepare data for clustering info
el01 <- CreateSeuratObject(counts = el01_toc, project = "Cont.Male")
el01 <- NormalizeData(el01, normalization.method = "LogNormalize", scale.factor = 10000)
el01 <- FindVariableFeatures(el01, selection.method = "vst", nfeatures = 2000)
el01 <- ScaleData(el01) # Scale data using all genes
el01 <- RunPCA(el01, npcs = 30, verbose = FALSE)
el01 <- RunUMAP(el01, reduction = "pca", dims = 1:10, min.dist = 0.5)
el01 <- FindNeighbors(el01, reduction = "pca", dims = 1:10)
el01 <- FindClusters(el01, resolution = 0.1)
DimPlot(el01, reduction="umap", label = T) + NoLegend() 

## Extract clustering informations
meta <- el01@meta.data
umap <- el01@reductions$umap@cell.embeddings

## Perform SoupX
el01_sc  <- setClusters(el01_sc, setNames(meta$seurat_clusters, rownames(meta)))
el01_sc  <- setDR(el01_sc, umap)
el01_sc  <- autoEstCont(el01_sc)
head(el01_sc$soupProfile[order(el01_sc$soupProfile$est, decreasing = T), ], n = 20)
adj.matrix  <- adjustCounts(el01_sc, roundToInt = T)

## Rerun the clustering procedure.
el01.soupx <- CreateSeuratObject(counts = adj.matrix, project = "Cont.Male", min.cells = 10, min.features = 200)
el01.soupx <- NormalizeData(el01.soupx, normalization.method = "LogNormalize", scale.factor = 10000)
el01.soupx <- FindVariableFeatures(el01.soupx, selection.method = "vst", nfeatures = 2000)
el01.soupx <- ScaleData(el01.soupx) # Scale data using all genes
el01.soupx <- RunPCA(el01.soupx, npcs = 30, verbose = FALSE)
el01.soupx <- RunUMAP(el01.soupx, reduction = "pca", dims = 1:10, min.dist = 0.5)
el01.soupx <- FindNeighbors(el01.soupx, reduction = "pca", dims = 1:10)
el01.soupx <- FindClusters(el01.soupx, resolution = 0.1)
DimPlot(el01.soupx, reduction="umap", label = T) + NoLegend() 
dim(el01.soupx) # 20580 9068
dim(el01) # 32285 9068

##### EL02 #####
## Prepare data for SoupX
el02_toc = Read10X(data.dir = "Rawdata/EL02_GEX/cellranger/EL02_GEX_cellranger_count_outs/filtered_feature_bc_matrix/")
el02_tod = Read10X(data.dir = "Rawdata/EL02_GEX/cellranger/EL02_GEX_cellranger_count_outs/raw_feature_bc_matrix/")
el02_sc = SoupChannel(el02_tod, el02_toc)
el02_sc

## Prepare data for clustering info
el02 <- CreateSeuratObject(counts = el02_toc, project = "Cont.Female")
el02 <- NormalizeData(el02, normalization.method = "LogNormalize", scale.factor = 10000)
el02 <- FindVariableFeatures(el02, selection.method = "vst", nfeatures = 2000)
el02 <- ScaleData(el02) # Scale data using all genes
el02 <- RunPCA(el02, npcs = 30, verbose = FALSE)
el02 <- RunUMAP(el02, reduction = "pca", dims = 1:10, min.dist = 0.5)
el02 <- FindNeighbors(el02, reduction = "pca", dims = 1:10)
el02 <- FindClusters(el02, resolution = 0.1)
DimPlot(el02, reduction="umap", label = T) + NoLegend() 

## Extract clustering informations
meta <- el02@meta.data
umap <- el02@reductions$umap@cell.embeddings

## Perform SoupX
el02_sc  <- setClusters(el02_sc, setNames(meta$seurat_clusters, rownames(meta)))
el02_sc  <- setDR(el02_sc, umap)
el02_sc  <- autoEstCont(el02_sc)
head(el02_sc$soupProfile[order(el02_sc$soupProfile$est, decreasing = T), ], n = 20)

adj.matrix  <- adjustCounts(el02_sc, roundToInt = T)

## Rerun the clustering procedure.
el02.soupx <- CreateSeuratObject(counts = adj.matrix, project = "Cont.Female", min.cells = 10, min.features = 200)
el02.soupx <- NormalizeData(el02.soupx, normalization.method = "LogNormalize", scale.factor = 10000)
el02.soupx <- FindVariableFeatures(el02.soupx, selection.method = "vst", nfeatures = 2000)
el02.soupx <- ScaleData(el02.soupx) # Scale data using all genes
el02.soupx <- RunPCA(el02.soupx, npcs = 30, verbose = FALSE)
el02.soupx <- RunUMAP(el02.soupx, reduction = "pca", dims = 1:10, min.dist = 0.5)
el02.soupx <- FindNeighbors(el02.soupx, reduction = "pca", dims = 1:10)
el02.soupx <- FindClusters(el02.soupx, resolution = 0.1)
DimPlot(el02.soupx, reduction="umap", label = T) + NoLegend() 
dim(el02) # 32285 7682
dim(el02.soupx) # 19780 7682

##### EL03 #####
## Prepare data for SoupX
el03_toc = Read10X(data.dir = "Rawdata/EL03_GEX/cellranger/EL03_GEX_cellranger_count_outs/filtered_feature_bc_matrix/")
el03_tod = Read10X(data.dir = "Rawdata/EL03_GEX/cellranger/EL03_GEX_cellranger_count_outs/raw_feature_bc_matrix/")
el03_sc = SoupChannel(el03_tod, el03_toc)
el03_sc

## Prepare data for clustering info
el03 <- CreateSeuratObject(counts = el03_toc, project = "Test.Male")
el03 <- NormalizeData(el03, normalization.method = "LogNormalize", scale.factor = 10000)
el03 <- FindVariableFeatures(el03, selection.method = "vst", nfeatures = 2000)
el03 <- ScaleData(el03) # Scale data using all genes
el03 <- RunPCA(el03, npcs = 30, verbose = FALSE)
el03 <- RunUMAP(el03, reduction = "pca", dims = 1:10, min.dist = 0.5)
el03 <- FindNeighbors(el03, reduction = "pca", dims = 1:10)
el03 <- FindClusters(el03, resolution = 0.1)
DimPlot(el03, reduction="umap", label = T) + NoLegend() 

## Extract clustering informations
meta <- el03@meta.data
umap <- el03@reductions$umap@cell.embeddings

## Perform SoupX
el03_sc  <- setClusters(el03_sc, setNames(meta$seurat_clusters, rownames(meta)))
el03_sc  <- setDR(el03_sc, umap)
el03_sc  <- autoEstCont(el03_sc)
head(el03_sc$soupProfile[order(el03_sc$soupProfile$est, decreasing = T), ], n = 20)

adj.matrix  <- adjustCounts(el03_sc, roundToInt = T)

## Rerun the clustering procedure.
el03.soupx <- CreateSeuratObject(counts = adj.matrix, project = "Test.Male", min.cells = 10, min.features = 200)
el03.soupx <- NormalizeData(el03.soupx, normalization.method = "LogNormalize", scale.factor = 10000)
el03.soupx <- FindVariableFeatures(el03.soupx, selection.method = "vst", nfeatures = 2000)
el03.soupx <- ScaleData(el03.soupx) # Scale data using all genes
el03.soupx <- RunPCA(el03.soupx, npcs = 30, verbose = FALSE)
el03.soupx <- RunUMAP(el03.soupx, reduction = "pca", dims = 1:10, min.dist = 0.5)
el03.soupx <- FindNeighbors(el03.soupx, reduction = "pca", dims = 1:10)
el03.soupx <- FindClusters(el03.soupx, resolution = 0.1)
DimPlot(el03.soupx, reduction="umap", label = T) + NoLegend() 
dim(el03)
dim(el03.soupx)


##### EL04 #####
## Prepare data for SoupX
el04_toc = Read10X(data.dir = "Rawdata/EL04_GEX/cellranger/EL04_GEX_cellranger_count_outs/filtered_feature_bc_matrix/")
el04_tod = Read10X(data.dir = "Rawdata/EL04_GEX/cellranger/EL04_GEX_cellranger_count_outs/raw_feature_bc_matrix/")
el04_sc = SoupChannel(el04_tod, el04_toc)
el04_sc

## Prepare data for clustering info
el04 <- CreateSeuratObject(counts = el04_toc, project = "Test.Female")
el04 <- NormalizeData(el04, normalization.method = "LogNormalize", scale.factor = 10000)
el04 <- FindVariableFeatures(el04, selection.method = "vst", nfeatures = 2000)
el04 <- ScaleData(el04) # Scale data using all genes
el04 <- RunPCA(el04, npcs = 30, verbose = FALSE)
el04 <- RunUMAP(el04, reduction = "pca", dims = 1:10, min.dist = 0.5)
el04 <- FindNeighbors(el04, reduction = "pca", dims = 1:10)
el04 <- FindClusters(el04, resolution = 0.1)
DimPlot(el04, reduction="umap", label = T) + NoLegend() 

## Extract clustering informations
meta <- el04@meta.data
umap <- el04@reductions$umap@cell.embeddings

## Perform SoupX
el04_sc  <- setClusters(el04_sc, setNames(meta$seurat_clusters, rownames(meta)))
el04_sc  <- setDR(el04_sc, umap)
el04_sc  <- autoEstCont(el04_sc)
head(el04_sc$soupProfile[order(el04_sc$soupProfile$est, decreasing = T), ], n = 20)

adj.matrix  <- adjustCounts(el04_sc, roundToInt = T)

## Rerun the clustering procedure.
el04.soupx <- CreateSeuratObject(counts = adj.matrix, project = "Test.Female", min.cells = 10, min.features = 200)
el04.soupx <- NormalizeData(el04.soupx, normalization.method = "LogNormalize", scale.factor = 10000)
el04.soupx <- FindVariableFeatures(el04.soupx, selection.method = "vst", nfeatures = 2000)
el04.soupx <- ScaleData(el04.soupx) # Scale data using all genes
el04.soupx <- RunPCA(el04.soupx, npcs = 30, verbose = FALSE)
el04.soupx <- RunUMAP(el04.soupx, reduction = "pca", dims = 1:10, min.dist = 0.5)
el04.soupx <- FindNeighbors(el04.soupx, reduction = "pca", dims = 1:10)
el04.soupx <- FindClusters(el04.soupx, resolution = 0.1)
DimPlot(el04.soupx, reduction="umap", label = T) + NoLegend() 
dim(el04)
dim(el04.soupx)

obj.list = list()
obj.list[["el01"]] = el01.soupx
obj.list[["el02"]] = el02.soupx
obj.list[["el03"]] = el03.soupx
obj.list[["el04"]] = el04.soupx

saveRDS(obj.list, "Seurat_Individual_SoupX_0622.RDS")

