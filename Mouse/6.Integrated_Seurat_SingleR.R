rm(list = ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
library(cowplot)
library(scRNAseq)
library(scuttle)
library(SingleR)
setwd("~/Dropbox/ANGPT2_AD_v2/")

combined = readRDS("Data/Seurat_Integrated_DecontX_AllBatch_Edit_0213.RDS")
combined.se = as.SingleCellExperiment(combined)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
## SingleR
### scRNAseq Reference
ziesel = ZeiselBrainData()
head(ziesel$level1class)
ziesel <- ziesel[,!is.na(ziesel$level1class)]
ziesel <- logNormCounts(ziesel)

tasic = TasicBrainData()
head(tasic$broad_type)
tasic <- tasic[,!is.na(tasic$broad_type)]
tasic <- logNormCounts(tasic)

## Ziesel
pred.combined <- SingleR(test = combined.se, ref = ziesel, de.method="wilcox",
                         labels = ziesel$level1class)
combined@meta.data$ziesel = pred.combined$pruned.labels

## Tasic
pred.combined <- SingleR(test = combined.se, ref = tasic, de.method="wilcox",
                         labels = tasic$broad_type)
combined@meta.data$tasic = pred.combined$pruned.labels

## TransferAnchor
allen_hipo <- readRDS('Resources/Allen_brain_expr_reference_seurat_object_onlyhipp_afternormpca.RDS')
anchors <- FindTransferAnchors(reference = allen_hipo, query = combined, 
                               normalization.method = "LogNormalize", dims = 1:25, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = allen_hipo$subclass_label, dims = 1:25)
combined <- AddMetaData(combined, metadata = predictions)
DimPlot(combined, reduction = "umap", group.by = "predicted.id", label = TRUE, pt.size = 0.3) #+ NoLegend()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

p1 = DimPlot(combined, label = TRUE, pt.size = .1) 
p2 = DimPlot(combined, group.by = "ziesel", label = TRUE, pt.size = .1) 
p3 = DimPlot(combined, group.by = "tasic", label = TRUE, pt.size = .1) 
p4 = DimPlot(combined, group.by = "predicted.id", label = TRUE, pt.size = .1) + guides(color = guide_legend(override.aes = list(size=4), ncol=2))

# p = plot_grid(p2, p3, p4, p1, ncol = 2)
p_1 = plot_grid(p2, p3)
p_2 = plot_grid(p4, p1, rel_widths = c(1.25, 1))
p = plot_grid(p_1, p_2, ncol = 1)
ggsave(p, file = "Figures/DimPlot_Seurat_Integrated_SingleR_Result.pdf", width = 14, height = 10)

saveRDS(combined, "Data/Seurat_Integrated_DecontX_Edit_SingleR_0215.RDS")
