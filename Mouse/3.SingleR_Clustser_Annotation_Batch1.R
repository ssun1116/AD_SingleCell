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

combined = readRDS("Data/data_combined_DecontX_Annotated_final.RDS")

comb1 = subset(combined, subset = orig.ident == "Control_F")
comb2 = subset(combined, subset = orig.ident == "Control_M")
comb3 = subset(combined, subset = orig.ident == "Test_F")
comb4 = subset(combined, subset = orig.ident == "Test_M")

sample.size = 7000
set.seed(111)

comb1 = subset(combined, subset = orig.ident == "Control_F")
orig.size1 = ncol(comb1) 
sampled.cells1 <- colnames(comb1)[sample(x = orig.size1, size = sample.size, replace = F)]
comb2 = subset(combined, subset = orig.ident == "Control_M")
orig.size2 = ncol(comb2) 
sampled.cells2 <- colnames(comb2)[sample(x = orig.size2, size = sample.size, replace = F)]
comb3 = subset(combined, subset = orig.ident == "Test_F")
orig.size3 = ncol(comb3) 
sampled.cells3 <- colnames(comb3)[sample(x = orig.size3, size = sample.size, replace = F)]
comb4 = subset(combined, subset = orig.ident == "Test_M")
orig.size4 = ncol(comb4) 
sampled.cells4 <- colnames(comb4)[sample(x = orig.size4, size = sample.size, replace = F)]

sampled.cells = c(sampled.cells1, sampled.cells2, sampled.cells3, sampled.cells4)
comb.sampled = subset(combined, cells = sampled.cells)

DimPlot(comb.sampled, group.by = "celltype", label = T)


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

p1 = DimPlot(combined, label = TRUE, pt.size = .1) + NoLegend()
p2 = DimPlot(combined, group.by = "ziesel", label = TRUE, pt.size = .1) 
p3 = DimPlot(combined, group.by = "tasic", label = TRUE, pt.size = .1) 
p4 = DimPlot(combined, group.by = "predicted.id", label = TRUE, pt.size = .1) + guides(color = guide_legend(override.aes = list(size=4), ncol=2))

# ggsave(p2, file = "Seurat_integrate_DecontX_Annotated_Ziesel.pdf", width = 7, height = 5)
# ggsave(p3, file = "Seurat_integrate_DecontX_Annotated_Tasic.pdf", width = 7, height = 5)
# ggsave(p4, file = "Seurat_integrate_DecontX_Annotated_Allen.pdf", width = 7, height = 5)

saveRDS(combined, "Data/data_combined_DecontX_Annotated.RDS")
combined = readRDS("Data/data_combined_DecontX_Annotated.RDS")

ggplot(combined@meta.data, aes(x=seurat_clusters, fill=ziesel)) + 
  geom_bar(position = "fill", color = "black") +
  theme_classic() +
  labs(x = '', y = 'Proportion of cell (%)') +
  theme(panel.grid.major = element_blank(),
        axis.text.x=element_text(size = 12, colour="black", hjust=0.5,vjust=0),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12, vjust = 0.2),
        # plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

ggplot(combined@meta.data, aes(x=seurat_clusters, fill=tasic)) + 
  geom_bar(position = "fill", color = "black") +
  theme_classic() +
  labs(x = '', y = 'Proportion of cell (%)') +
  theme(panel.grid.major = element_blank(),
        axis.text.x=element_text(size = 12, colour="black", hjust=0.5,vjust=0),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12, vjust = 0.2),
        # plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

ggplot(combined@meta.data, aes(x=seurat_clusters, fill=predicted.id)) + 
  geom_bar(position = "fill", color = "black") +
  theme_classic() +
  labs(x = '', y = 'Proportion of cell (%)') +
  theme(panel.grid.major = element_blank(),
        axis.text.x=element_text(size = 12, colour="black", hjust=0.5,vjust=0),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12, vjust = 0.2),
        # plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

# 
# ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# ### Final Annotation
# combined = readRDS("Data/data_combined_DecontX_Annotated_0511.RDS")
# DimPlot(combined, label = T)
# 
# annot.table = data.frame(seurat_clusters = c(0:17),
#                          celltype = c("Oligo", "CR", "Ex", "Ex", "In", "Ex", "Ex.DG", "Astro", "Micro",
#                                       "In", "EC.MC", "Ex", "Ex.CA1", "OPC", "In", "VLMC", "Oligo", "CR"))
# meta = as.data.frame(combined@meta.data)
# meta$barcode = rownames(meta)
# meta2 = merge(meta, annot.table, by = "seurat_clusters")
# rownames(meta2) = meta2$barcode
# combined <- AddMetaData(combined, meta2)
# 
# DimPlot(combined, group.by = "celltype", label = T)
# 
# ### Need to further process the EC group
# combined$celltype = ifelse(!is.na(combined$tasic), ifelse(combined$tasic == "Astrocyte", "Astro", combined$celltype), combined$celltype)
# combined$celltype = ifelse(!is.na(combined$predicted.id), ifelse(combined$predicted.id == "VLMC", "VLMC", combined$celltype), combined$celltype)
# combined$celltype = ifelse(combined$celltype == "EC.MC" & !is.na(combined$predicted.id), 
#                           ifelse(combined$predicted.id %in% c("CR", "Oligo"), "CR", combined$celltype), combined$celltype)
# 
# DefaultAssay(combined) = "RNA"
# combined = SetIdent(combined, value = "celltype")
# DimPlot(combined, group.by = "celltype", label = T)
# 
# combined$cond = ifelse(combined$orig.ident %in% c("Control_M", "Control_F"), "Control", "Test")
# combined$sex = ifelse(combined$orig.ident %in% c("Control_M", "Test_M"), "Male", "Female")
# 
# saveRDS(combined, "Data/data_combined_DecontX_Annotated_final.RDS")
# ggsave(last_plot(), file = "Seurat_integrate_DecontX_Annotated_Final.pdf", width = 6, height = 5)
# 
# p1 = DimPlot(combined, group.by = "celltype", label = T)
# p2 = DimPlot(combined, group.by = "orig.ident")
# library(cowplot)
# plot_grid(p1, p2, ncol = 2)
# 
# 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# basic plot of clusters by replicate
combined_F = subset(comb.sampled, subset = orig.ident %in% c("Control_F", "Test_F"))
combined_M = subset(comb.sampled, subset = orig.ident %in% c("Control_M", "Test_M"))

ggplot(comb.sampled@meta.data, aes(x=celltype, fill=orig.ident)) +
  geom_bar(position = "fill", color = "black") +
  theme_classic() +
  labs(x = '', y = 'Proportion of cell (%)') +
  theme(panel.grid.major = element_blank(),
        axis.text.x=element_text(size = 12, colour="black", hjust=1,vjust=1, angle = 45),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12, vjust = 0.2),
        # plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

ggplot(combined_F@meta.data, aes(x=celltype, fill=orig.ident)) +
  geom_bar(position = "fill", color = "black") +
  theme_classic() +
  labs(x = '', y = 'Proportion of cell (%)') +
  theme(panel.grid.major = element_blank(),
        axis.text.x=element_text(size = 12, colour="black", hjust=1,vjust=1, angle = 45),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12, vjust = 0.2),
        # plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

p1.2 = ggplot(combined_M@meta.data, aes(x=celltype, fill=orig.ident)) +
  geom_bar(position = "fill", color = "black") +
  theme_classic() +
  labs(x = '', y = 'Proportion of cell (%)') +
  theme(panel.grid.major = element_blank(),
        axis.text.x=element_text(size = 12, colour="black", hjust=1,vjust=1, angle = 45),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12, vjust = 0.2),
        # plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

plot_grid(p1.1, p1.2)

p2 = ggplot(combined@meta.data, aes(x=celltype, fill=cond)) +
  geom_bar(position = "fill", color = "black") +
  theme_classic() +
  labs(x = '', y = 'Proportion of cell (%)') +
  theme(panel.grid.major = element_blank(),
        axis.text.x=element_text(size = 12, colour="black", hjust=1,vjust=1, angle = 45),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12, vjust = 0.2),
        # plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  geom_hline(yintercept = 0.5, colour = "black", size = 1)


plot_grid(p1, p2)

FeaturePlot(combined, features = c("Nrgn", "Gad2", "Aqp4", "Mbp", "Csf1r", "Vcan", "Flt1", "Tagln", "Atp13a5", "Cemip")) +
  patchwork::plot_layout(ncol = 5, nrow = 2)




