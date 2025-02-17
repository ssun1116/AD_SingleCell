rm(list = ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(writexl)

combined = readRDS("Data/Seurat_Integrated_DecontX_Edit_SingleR.RDS")
DimPlot(combined, label = T) # seurat_clusters resolution 0.7

DefaultAssay(combined) = "RNA"
d = FindAllMarkers(combined, logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE)
write_xlsx(d, path = "Tables/FindAllMarkers_Seurat_Integated_by_cluster_res.0.7.xlsx", 
           col_names = TRUE)

combined <- BuildClusterTree(object = combined)
PlotClusterTree(object = combined)

######
DefaultAssay(combined) = "RNA"
p = VlnPlot(combined, features = c("Slc17a7", "Gad1"), pt.size = 0) &
  theme(axis.text.x = element_text(angle= 0, hjust = 0.5),
        axis.title.x = element_blank())
ggsave(p, filename = "Figures/VlnPlot_Seurat_Integrated_Slc17a7.Gad1.pdf",
        width = 12, height= 4)

FeaturePlot(combined, features = c("Eya2", "Foxd1"))
VlnPlot(combined, features = c("Eya2", "Foxd1")) &
  theme(axis.text.x = element_text(angle= 0, hjust = 0.5),
        axis.title.x = element_blank())

## Filter out VLMC cells
combined2 = subset(combined, !(combined$seurat_clusters == "12" & (Eya2 > 1 | Foxd1 > 1)))

VlnPlot(combined2, features = c("Eya2", "Foxd1"))
DimPlot(combined2, label = T)
FeaturePlot(combined2, features = "Dcn")

## Hippocampus markers
FeaturePlot(combined2, features = c("Plekhg1", "Lefty1", "Gpc3", "Dlk1", "Wfs1"), ncol = 5) ## CA1
FeaturePlot(combined2, features = c("Ramp3",  "Map3k15")) ## CA2
FeaturePlot(combined, features = c("Iyd", "Fam197a", "Pvrl3")) ## CA3
FeaturePlot(combined, features = c("Trp73")) ## CR
FeaturePlot(combined2, features = c("Prox1", "Crlf1")) # DG

### BICCN markers
FeaturePlot(combined, features = c("Sv2b", "Arpp21", "Gad1", "Dlx6os1", "Gli3", "Slc39a12", "Gjc3", "Pdgfra", "Stk32a", "Inpp5d", "Adgrl4", "Slco1a4", "Eya2", "Foxd1", "Acta2", "Abcc9", "Tbx3os1"), ncol = 5)

ggsave(last_plot(), filename = "Figures/FeaturePlot_Seurat_Integrated_BICCN_Marker.pdf", width = 15, height = 10)





### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Final Annotation
combined2$celltype <- NULL # Combined2: VLMC excluded data
annot.table = data.frame(seurat_clusters = c(0:23),
                         celltype = c("Oligo", "CR", "Inh-1", "Ext-1", "Inh-2", "Ext-2", "Oligo", "Ext-3", "Ext-DG", "Ext-2", "Astro", "Micro", "Endo", "Ext-4", "Ext-5", "Astro", "OPC", "Micro", "Inh-1", "Ext-6", "Low-quality", "Micro", "Ext-DG", "Ext-DG"))
meta = as.data.frame(combined2@meta.data)
meta$barcode = rownames(meta)
meta2 = merge(meta, annot.table, by = "seurat_clusters")
rownames(meta2) = meta2$barcode
combined2 <- AddMetaData(combined2, meta2)

DimPlot(combined2, group.by = "celltype", label = T)

combined2$celltype = factor(combined2$celltype, levels = c("Ext-1", "Ext-2", "Ext-3", "Ext-4", "Ext-5","Ext-6", "Ext-DG", "Inh-1", "Inh-2", "Astro",  "OPC", "Oligo","Micro", "CR", "Endo", "Low-quality"))
combined2 = SetIdent(combined2, value = "celltype")

my_cols <- c("Ext-1" = "#7CAE00", "Ext-2" = "#57A101", "Ext-3" = "#319601", "Ext-4" = "#0CB702", "Ext-5" = "#00Be67", "Ext-6" = "#00C19A", "Ext-7" = "#17fc09", "Ext-DG" = "#00BFC4", "Inh-1" = "#00B8E7", "Inh-2" = "#00A9FF", "Astro" = "#F8766D", "Micro" = "#C77CFF", "Oligo" = "#ED68ED", "OPC" = "#FF68A1", "CR" = "#E68613", "Endo" = "#CD9600", "Low-quality" = "grey70")

# combined3 = subset(combined2, subset = celltype != "Low-quality")

saveRDS(combined2, "Data/Seurat_Integrated_DecontX_Edit_DraftAnnot_0219.RDS")

DimPlot(combined2, cols = my_cols, label=TRUE) 

ggsave(last_plot(), file = "Figures/DimPlot_Seurat_Integrated_DecontX_Edit_DraftAnnot_0219.pdf", width = 6, height = 5)

VlnPlot(combined2, features = c("Slc17a7", "Gad1"), ncol = 2, pt.size = 0)

d = FindAllMarkers(combined2, logfc.threshold = 0.25, min.pct = 0.1, only.pos = TRUE)
write_xlsx(d, path = "Tables/FindAllMarkers_Seurat_Integated_by_celltype.xlsx", 
           col_names = TRUE)

d_top = d %>% filter(cluster %in% c("Ext-1", "Ext-2", "Ext-3", "Ext-4", "Ext-DG", "Inh-1", "Inh-2", "CR", "Low-quality")) %>% group_by(cluster) %>%top_n(10, abs(avg_log2FC)) %>% pull(gene) %>% unique()

DotPlot(combined2, features = d_top) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), axis.title.y = element_blank()) + NoLegend()



# 
# ### Ziesel markers
# FeaturePlot(combined, features = c("Thy1", "Gad1", "Tbr1", "Spink8", "Mbp", "Aldoc", "Aif1", "Cldn5", "Acta2"), ncol = 4)
#
# FeaturePlot(combined, features = c("Lamp5", "Sncg", "Vip", "Sst", "Pvalb", "Meis2"))
# FeaturePlot(combined, features = c("Reln"))
# 
# ### Dotplot - unclarified clusters
# features1 = c("Rarb", "Rgs9", "Phactr1", "Pde10a", "Ryr3", "Meis2", "Pde7b")
# features2 = c("Cmss1", "Ttr", "Gapdh", "Calm1", "Nrgn", "Jarid2", "Camk2n1")
# features3 = c("Epha6", "Tafa1", "Shisa6", "Galnt16", "Iqgap2", "Cntnap5c", "Galnt17")
# features4 = c("Ntng1", "Synpo2", "Kcnc2", "Vav3", "Inpp4b", "Nell1", "Hdac9")
# features5 = c("Epha7", "Glis3", "Pcsk5", "Tafa2", "Ppfia2", "Rfx3")
# 
# p_ex = DotPlot(combined, features = c(features3, features5)) + 
#   theme(axis.title.x = element_blank())
# 
# p_in = DotPlot(combined, features = c(features1, features4)) + 
#   theme(axis.title.x = element_blank())
# 
# p_cr = DotPlot(combined, features = c(features2)) + 
#   theme(axis.title.x = element_blank())
# 
# ggsave(p_ex, filename = "Figures/DotPlot_Seurat_Integrate_Excitatory_8_14_marker.pdf", width = 12, height = 4)
# 
# ggsave(p_in, filename = "Figures/DotPlot_Seurat_Integrate_Inhibitory_4_13_marker.pdf", width = 12, height = 4)
# 
# ggsave(p_cr, filename = "Figures/DotPlot_Seurat_Integrate_CR_1_marker.pdf", width = 8, height = 4)
# 