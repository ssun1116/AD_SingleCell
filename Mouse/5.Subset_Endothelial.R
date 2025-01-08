library(Seurat)
library(dplyr)
library(RColorBrewer)

combined = readRDS("Data/data_combined_DecontX_Annotated_final.RDS")
endo = subset(combined, subset = celltype =="EC.MC")
table(endo$celltype)
DimPlot(endo, label = T)

endo.list = list()
endo.list[["EL01"]] = subset(endo, subset = orig.ident == "Control_M")
endo.list[["EL02"]] = subset(endo, subset = orig.ident == "Control_F")
endo.list[["EL03"]] = subset(endo, subset = orig.ident == "Test_M")
endo.list[["EL04"]] = subset(endo, subset = orig.ident == "Test_F")

## Harmony
merged2 = Reduce(merge, endo.list)
merged2 <- NormalizeData(merged2, normalization.method = "LogNormalize", scale.factor = 10000)
merged2 <- FindVariableFeatures(merged2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(merged2) # Scale data using all genes
merged2 <- ScaleData(merged2, features = all.genes)
merged2 <- RunPCA(merged2, features = VariableFeatures(object = merged2))
merged2 <- merged2 %>% RunHarmony("orig.ident", plot_convergence = TRUE) %>% 
  RunUMAP(reduction = "harmony", dims = 1:5, min.dist = 0.3) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:5) %>% 
  FindClusters(resolution = 0.4) 

DimPlot(object = merged2, reduction = "umap", label = T)

merged2$celltype2 = ifelse(merged2$seurat_clusters == 0, "Pericyte",
                           ifelse(merged2$seurat_clusters == 1, "Endothelial",
                                  ifelse(merged2$seurat_clusters == 2, "Mural", "Immune")))

saveRDS(merged2, "Data/data_combined_DecontX_Endothelial.RDS")

d = FindAllMarkers(merged2, only.pos = TRUE)

features = c("Vegfc", "Dkk2", "Arl15", "Pecam1", "Flt1", "Cldn5", "Atp10a", "Syne1", "Abcb1", "Npipb5", "Cmtm8", "Slc2a1",
             "Tshz2", "Adgrg6", "Aff3", "Grm8", "Pdgfrb", "Cers6", "Dab1", "Shisa6", "Itm2a", "Ca4", "Spock2", "Slc38a5",
             "Acta2", "Mcam", "Hbb", "Myocd", "Myh11", "Kcnab1", "Ntn4", "Gfpt2", "Col12a1", "Cemip", "Svil", "Rnf220",
             "Anxa2", "Slc13a3", "Trpm3", "Tp63", "Kcnma1", "Slit2", "Eys")

DoHeatmap(endo, features = features)
FeaturePlot(merged2, features =c("Vegfc", "Arl15", "Mfsd2a", "Slc7a5", "Tshz2", "Adgrg6"))
FeaturePlot(merged2, features = c("Acta2", "Myh11", "Myocd", "Cd74"))

deg.list = list()
for (cell in unique(merged2$celltype2)){
  merged2.sub = subset(merged2, subset = celltype2 == cell)
  deg = FindMarkers(merged2.sub, group.by = "cond", ident.1 = "Test", ident.2 = "Control", logfc.threshold = -Inf)
  deg$gene = rownames(deg)
  deg$change = ifelse(deg$p_val_adj < 0.05 & abs(deg$avg_log2FC) > 0.25, 
                      ifelse(deg$avg_log2FC > 0, "Up", "Down"), "None")
  deg$cluster = cell
  deg.list[[cell]] = deg
}
deg_by_clust <- bind_rows(deg.list)
deg_by_clust$gene <- gsub("\\..*", "", rownames(deg_by_clust))

deg_by_clust2 = deg_by_clust %>% filter(change != "None")
ggplot(deg_by_clust2, aes(x = cluster, fill = change)) +  # Plot with values on top
  geom_bar(alpha = 0.9) +
  scale_fill_manual(values = c("#3182bd", "#E64B35")) +
  theme_classic() +
  labs(x = '', y = 'Number of DEGs') +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x=element_text(size = 12, colour="black", hjust=0.5,vjust=0),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.y = element_text(size = 12, vjust = 0.2),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) + coord_flip()

degs = FindMarkers(merged2, group.by = "cond", ident.1 = "Test", ident.2 = "Control", logfc.threshold = -Inf)



