rm(list = ls())
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(SCP)

combined = readRDS("Data/Seurat_Integrated_DecontX_Edit_DraftAnnot_0219.RDS")

# d1 = srt_to_adata(combined, assay_X = "RNA", slot_layers = "data") # run in imac.
# d1$write_h5ad("Data/Seurat_Integrated_DecontX_Edit_DraftAnnot_0219.h5ad")

d1 = srt_to_adata(combined, assay_X = "RNA", slot_layers = "counts") # run in imac.
d1$write_h5ad("Data/Seurat_Integrated_DecontX_Edit_DraftAnnot_counts_0402.h5ad")


