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

combined = readRDS("Data/data_combined_SoupX_Doublet_Sampled_0626.RDS")
DimPlot(combined, group.by = "orig.ident")

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


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

p1 = DimPlot(combined, label = TRUE, pt.size = .1) 
p2 = DimPlot(combined, group.by = "ziesel", label = TRUE, pt.size = .1)
p3 = DimPlot(combined, group.by = "tasic", label = TRUE, pt.size = .1) 
p4 = DimPlot(combined, group.by = "predicted.id", label = TRUE, pt.size = .1) + 
  guides(color = guide_legend(override.aes = list(size=4), ncol=2))

plot_grid(p1, p2, p3, p4, ncol = 2)

# ggsave(p2, file = "Seurat_integrate_SoupX_Annotated_Ziesel.pdf", width = 7, height = 5)
# ggsave(p3, file = "Seurat_integrate_SoupX_Annotated_Tasic.pdf", width = 7, height = 5)
# ggsave(p4, file = "Seurat_integrate_SoupX_Annotated_Allen.pdf", width = 7, height = 5)

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


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Final Annotation

annot.table = data.frame(seurat_clusters = c(0:18),
                         celltype = c("Oligo", "CR", "Ex1", "In", "Ex1", "Ex.DG", "Ex3", "Astro", "Ex1", "Micro",
                                      "Ex2", "Ex1", "OPC", "EC.MC", "CR", "Oligo", "Oligo", "VLMC", "Ex1"))
meta = as.data.frame(combined@meta.data)
meta$barcode = rownames(meta)
meta2 = merge(meta, annot.table, by = "seurat_clusters")
rownames(meta2) = meta2$barcode
combined <- AddMetaData(combined, meta2)

DimPlot(combined, group.by = "celltype", label = T)

combined$cond = ifelse(combined$orig.ident %in% c("Cont.Male", "Cont.Female"), "Control", "Test")
combined$sex = ifelse(combined$orig.ident %in% c("Cont.Male", "Test.Male"), "Male", "Female")

FeaturePlot(combined, features = c("Nrgn", "Gad2", "Aqp4", "Mbp", "Csf1r", "Vcan", "Flt1", "Tagln", "Atp13a5", "Cemip")) +
  patchwork::plot_layout(ncol = 5, nrow = 2)

combined = SetIdent(combined, value = "celltype")

combined <- FindSubCluster(combined, "Ex1", graph.name = "integrated_snn", subcluster.name = "sub.cluster", resolution= 0.5)
DimPlot(combined, group.by = "sub.cluster")

combined$celltype = ifelse(combined$sub.cluster == "Ex1_5", "In", combined$celltype)
DimPlot(combined, group.by = "celltype", label = T)

saveRDS(combined, "Data/data_combined_SoupX_Doublet_Annotated_0627.RDS")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# basic plot of clusters by replicate

p1 = ggplot(combined@meta.data, aes(x=celltype, fill=orig.ident)) +
  geom_bar(position = "fill", color = "black") +
  theme_classic() +
  labs(x = '', y = 'Proportion of cell (%)') +
  theme(panel.grid.major = element_blank(),
        axis.text.x=element_text(size = 12, colour="black", hjust=1,vjust=1, angle = 45),
        axis.text.y = element_text(size = 10, colour = "black"),
        axis.title.y = element_text(size = 12, vjust = 0.2),
        # plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

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


# advanced plot of clusters by sex
combined_F = subset(combined, subset = sex == "Female")
combined_M = subset(combined, subset = sex != "Female")

p1.1 = ggplot(combined_F@meta.data, aes(x=celltype, fill=orig.ident)) +
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

p1.2 = ggplot(combined_M@meta.data, aes(x=celltype, fill=orig.ident)) +
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

plot_grid(p1.1, p1.2, labels = c("Female", "Male"))


