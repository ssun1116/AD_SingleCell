rm(list = ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(writexl)

combined = readRDS("Data/Seurat_Integrated_DecontX_Edit_DraftAnnot_0219.RDS")
DimPlot(combined, label = T)
DimPlot(combined, group = "seurat_clusters", label = T)

                 
             
mapmycells = read.csv("MapMyCells_10xWholeMouseBrain_HierarchicalMapping_Result.csv")
mapmycells2 = mapmycells[,c("cell_id", "class_name")]
rownames(mapmycells2) = mapmycells2$cell_id
mapmycells2$cell_id <- NULL

meta = as.data.frame(combined@meta.data)
meta2 = cbind.data.frame(meta, mapmycells2)

combined@meta.data = meta2

DimPlot(combined, group = "class_name", label = T)

ggplot(combined@meta.data, aes(x=seurat_clusters, fill=class_name)) + 
  geom_bar(position = "fill", color = "black") +
  theme_classic() +
  labs(x = '', y = 'Proportion of cell (%)') +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x=element_text(size = 12, colour="black", angle = 45, hjust=1,vjust=1),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.y = element_text(size = 12, vjust = 0.2),
    # plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))



annot.table = data.frame(seurat_clusters = c(0:23),
                         final_annot = c("Oligo", "ChP", "Inh (CGE/MGE)", "Ext (Cortex)", "Inh (LGE)", "Ext (Cortex)", "Oligo", "Ext (Cortex)", 
                                      "Ext (Hippo)", "Ext (Cortex)", "Astro", "Micro", "Endo", "Ext (Thalamus)", "Ext (Hippo)", "Astro", "OPC",
                                      "Micro", "Inh (CGE/MGE)", "Ext (Thalamus)", "Ependymal", "Micro", "Ext (Hippo)", "Ext (Hippo)"))
meta = as.data.frame(combined@meta.data)
meta$barcode = rownames(meta)
meta2 = merge(meta, annot.table, by = "seurat_clusters")
rownames(meta2) = meta2$barcode
combined <- AddMetaData(combined, meta2)

DimPlot(combined, group.by = "final_annot", label = T)

# my_cols <- c("Ext-1" = "#7CAE00", "Ext-2" = "#57A101", "Ext-3" = "#319601", 
#              "Ext-4" = "#0CB702", "Ext-5" = "#00Be67", "Ext-6" = "#00C19A", 
#              "Ext-7" = "#17fc09", "Ext-DG" = "#00BFC4", "Inh-1" = "#00B8E7", 
#              "Inh-2" = "#00A9FF", "Astro" = "#F8766D", "Micro" = "#C77CFF", 
#              "Oligo" = "#ED68ED", "OPC" = "#FF68A1", "CR" = "#E68613", 
#              "Endo" = "#CD9600", "Low-quality" = "grey70")

my_cols <- c("Ext (Cortex)" = "#7CAE00", "Ext (Hippo)" = "#00Be67", "Ext (Thalamus)" = "#319601", 
             "Inh (CGE/MGE)" = "#00BFC4", "Inh (LGE)" = "#00A9FF", 
             "Astro" = "#F8766D", "Micro" = "#C77CFF", "Oligo" = "#ED68ED", "OPC" = "#FF68A1", 
             "ChP" = "#E68613", "Endo" = "#CD9600", "Ependymal" = "#FEC20C")

# DimPlot(combined, group.by = "final_annot", cols = my_cols, label = T)

combined$final_annot = factor(combined$final_annot, 
                              levels = c("Ext (Cortex)", "Ext (Hippo)", "Ext (Thalamus)", 
                                         "Inh (CGE/MGE)", "Inh (LGE)",
                                         "Astro",  "OPC", "Oligo","Micro", "ChP", "Endo", "Ependymal"))

combined = SetIdent(combined, value = "final_annot")


## 240827 update.

combined = readRDS("Data/Seurat_Integrated_DecontX_Edit_DraftAnnot_0519.RDS")
combined$final_annot = as.character(combined$final_annot)
combined$final_annot = ifelse(combined$final_annot == "Endo", "Vascular", combined$final_annot)

my_cols <- c("Ext (Cortex)" = "#7CAE00", "Ext (Hippo)" = "#00Be67", "Ext (Thalamus)" = "#319601", 
             "Inh (CGE/MGE)" = "#00BFC4", "Inh (LGE)" = "#00A9FF", 
             "Astro" = "#F8766D", "Micro" = "#C77CFF", "Oligo" = "#ED68ED", "OPC" = "#FF68A1", 
             "ChP" = "#E68613", "Vascular" = "#CD9600", "Ependymal" = "#FEC20C")

# DimPlot(combined, group.by = "final_annot", cols = my_cols, label = T)

combined$final_annot = factor(combined$final_annot, 
                              levels = c("Ext (Cortex)", "Ext (Hippo)", "Ext (Thalamus)", 
                                         "Inh (CGE/MGE)", "Inh (LGE)",
                                         "Astro",  "OPC", "Oligo","Micro", "ChP", "Vascular", "Ependymal"))

combined = SetIdent(combined, value = "final_annot")

p = DimPlot(combined, group.by = "final_annot", cols = my_cols, label = T)

# saveRDS(combined, "Data/Seurat_Integrated_DecontX_Edit_DraftAnnot_0519.RDS")
saveRDS(combined, "Data/Seurat_Integrated_DecontX_Edit_DraftAnnot_0827.RDS")

ggsave("Figures/ANGPT2_AD/DimPlot_Seurat_Integrated_DecontX_Edit_DraftAnnot_0827.pdf", p, width = 6, height = 4.5)

### Visualization



d = FindAllMarkers(combined, logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE)
write_xlsx(d, path = "Tables/FindAllMarkers_Seurat_Integated_by_final_annot_0519.xlsx", col_names = TRUE)

p = FeaturePlot(combined, features = c("Satb2", "Tbr1", "Shisa6", "Epha6", "Slc17a6", "Tcf7l2", "Adarb2", "Lhx6", "Meis2", "Ebf1", 
                                   "Ttr", "Dkk3", "Foxj1", "Pifo", "Pdgfra", "Gli3", "Gjc3", "Inpp5d", "Adgrl4"), ncol = 5, raster = TRUE)

ggsave("Figures/ANGPT2_AD/FeaturePlot_Seurat_Integrated_DecontX_Edit_DraftAnnot_raster_1118.pdf", p, width = 16, height = 12)

p = FeaturePlot(combined, features = c("Satb2", "Tbr1", "Shisa6", "Epha6", "Slc17a6", "Tcf7l2", "Adarb2", "Lhx6", "Meis2", "Ebf1", 
                                       "Ttr", "Dkk3", "Foxj1", "Pifo", "Pdgfra", "Gli3", "Gjc3", "Inpp5d", "Adgrl4"), ncol = 5, raster = FALSE)

ggsave("Figures/ANGPT2_AD/FeaturePlot_Seurat_Integrated_DecontX_Edit_DraftAnnot_1118.pdf", p, width = 16, height = 12)

