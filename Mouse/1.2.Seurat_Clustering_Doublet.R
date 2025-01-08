rm(list = ls())
library(celda)
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(DoubletFinder)
library(SingleCellExperiment)
setwd("~/Dropbox/ANGPT2_AD/")

obj.list = readRDS("Data/Seurat_Individual_SoupX_0622.RDS")
seu_el01 = obj.list[["el01"]]
seu_el02 = obj.list[["el02"]]
seu_el03 = obj.list[["el03"]]
seu_el04 = obj.list[["el04"]]

##### EL01 #####

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
## ground-truth: antibody-barcode sample multiplexing (Cell Hashing) and SNP deconvolution (Demuxlet)
sweep.res.list_el01 <- paramSweep_v3(seu_el01, PCs = 1:10, sct = FALSE)
sweep.stats_el01 <- summarizeSweep(sweep.res.list_el01, GT = FALSE)
bcmvn_el01 <- find.pK(sweep.stats_el01)
bcmvn_el01 = as.data.frame(bcmvn_el01)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seu_el01@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_el01@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(seu_el01@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu_el01 <- doubletFinder_v3(seu_el01, PCs = 1:10, pN = 0.25, pK = as.numeric(bcmvn_el01[which.max(bcmvn_el01$BCmetric),]$pK), 
                             nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu_el01 <- doubletFinder_v3(seu_el01, PCs = 1:10, pN = 0.25, pK = as.numeric(bcmvn_el01[which.max(bcmvn_el01$BCmetric),]$pK), 
                             nExp = nExp_poi.adj, reuse.pANN = colnames(seu_el01@meta.data)[6], sct = FALSE)


##### el02 #####

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_el02 <- paramSweep_v3(seu_el02, PCs = 1:10, sct = FALSE)
sweep.stats_el02 <- summarizeSweep(sweep.res.list_el02, GT = FALSE)
bcmvn_el02 <- find.pK(sweep.stats_el02)
bcmvn_el02 = as.data.frame(bcmvn_el02)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seu_el02@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_el02@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(seu_el02@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu_el02 <- doubletFinder_v3(seu_el02, PCs = 1:10, pN = 0.25, pK = as.numeric(bcmvn_el02[which.max(bcmvn_el02$BCmetric),]$pK), 
                             nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu_el02 <- doubletFinder_v3(seu_el02, PCs = 1:10, pN = 0.25, pK = as.numeric(bcmvn_el02[which.max(bcmvn_el02$BCmetric),]$pK), 
                             nExp = nExp_poi.adj, reuse.pANN = colnames(seu_el02@meta.data)[6], sct = FALSE)


##### el03 #####

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_el03 <- paramSweep_v3(seu_el03, PCs = 1:10, sct = FALSE)
sweep.stats_el03 <- summarizeSweep(sweep.res.list_el03, GT = FALSE)
bcmvn_el03 <- find.pK(sweep.stats_el03)
bcmvn_el03 = as.data.frame(bcmvn_el03)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seu_el03@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_el03@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(seu_el03@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu_el03 <- doubletFinder_v3(seu_el03, PCs = 1:10, pN = 0.25, pK = as.numeric(bcmvn_el03[which.max(bcmvn_el03$BCmetric),]$pK), 
                             nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu_el03 <- doubletFinder_v3(seu_el03, PCs = 1:10, pN = 0.25, pK = as.numeric(bcmvn_el03[which.max(bcmvn_el03$BCmetric),]$pK), 
                             nExp = nExp_poi.adj, reuse.pANN = colnames(seu_el03@meta.data)[6], sct = FALSE)


##### el04 #####

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_el04 <- paramSweep_v3(seu_el04, PCs = 1:10, sct = FALSE)
sweep.stats_el04 <- summarizeSweep(sweep.res.list_el04, GT = FALSE)
bcmvn_el04 <- find.pK(sweep.stats_el04)
bcmvn_el04 = as.data.frame(bcmvn_el04)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seu_el04@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_el04@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(seu_el04@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu_el04 <- doubletFinder_v3(seu_el04, PCs = 1:10, pN = 0.25, pK = as.numeric(bcmvn_el04[which.max(bcmvn_el04$BCmetric),]$pK), 
                             nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu_el04 <- doubletFinder_v3(seu_el04, PCs = 1:10, pN = 0.25, pK = as.numeric(bcmvn_el04[which.max(bcmvn_el04$BCmetric),]$pK), 
                             nExp = nExp_poi.adj, reuse.pANN = colnames(seu_el04@meta.data)[6], sct = FALSE)

obj.list = list()
obj.list[["el01"]] = seu_el01
obj.list[["el02"]] = seu_el02
obj.list[["el03"]] = seu_el03
obj.list[["el04"]] = seu_el04

saveRDS(obj.list, "Data/Seurat_Individual_SoupX_Doublet_0626.RDS")


