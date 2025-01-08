rm(list = ls())
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

combined = readRDS("Data/data_combined_DecontX_Annotated_final.RDS")
DimPlot(combined, group.by = "celltype", label = T)
head(combined@meta.data)
table(combined$orig.ident)
table(combined$cond)

colnames(combined@assays$RNA)

deg.list = list()
for (cell in unique(combined$celltype)){
  deg = FindMarkers(combined, subset.ident = cell,
                    group.by = "cond", ident.1 = "Test", ident.2 = "Control", 
                    test.use = "LR", logfc.threshold = -Inf)
  deg$gene = rownames(deg)
  deg$change = ifelse(deg$p_val_adj < 0.05, ifelse(deg$avg_log2FC > 0, "Up", "Down"), "None")
  deg$cluster = cell
  deg.list[[cell]] = deg
}
deg_by_clust <- bind_rows(deg.list)
deg_by_clust$gene <- gsub("\\..*", "", rownames(deg_by_clust))

saveRDS(deg_by_clust, file = "Results/deg_by_clust_table_DecontX_final_DEG.All.RDS")

deg.female = list()
for (cell in unique(combined$celltype)){
  deg = FindMarkers(combined, subset.ident = cell,
                    group.by = "cond.sex", ident.1 = "Test.Female", ident.2 = "Control.Male", 
                    test.use = "LR", logfc.threshold = -Inf)
  deg$gene = rownames(deg)
  deg$change = ifelse(deg$p_val_adj < 0.05, ifelse(deg$avg_log2FC > 0, "Up", "Down"), "None")
  deg$cluster = cell
  deg.female[[cell]] = deg
}
deg_female <- bind_rows(deg.female)
deg_female$gene <- gsub("\\..*", "", rownames(deg_female))

deg.male = list()
for (cell in unique(combined$celltype)){
  deg = FindMarkers(combined, subset.ident = cell,
                    group.by = "cond.sex", ident.1 = "Test.Female", ident.2 = "Control.Male", 
                    test.use = "LR", logfc.threshold = -Inf)
  deg$gene = rownames(deg)
  deg$change = ifelse(deg$p_val_adj < 0.05, ifelse(deg$avg_log2FC > 0, "Up", "Down"), "None")
  deg$cluster = cell
  deg.male[[cell]] = deg
}
deg_male <- bind_rows(deg.male)
deg_male$gene <- gsub("\\..*", "", rownames(deg_male))

deg.sex = list()
for (cell in unique(combined$celltype)){
  deg = FindMarkers(combined, subset.ident = cell,
                    group.by = "cond", ident.1 = "Test", ident.2 = "Control", 
                    test.use = "LR", logfc.threshold = -Inf, vars.to.regress = "sex")
  deg$gene = rownames(deg)
  deg$change = ifelse(deg$p_val_adj < 0.05, ifelse(deg$avg_log2FC > 0, "Up", "Down"), "None")
  deg$cluster = cell
  deg.sex[[cell]] = deg
}
deg_sex <- bind_rows(deg.sex)
deg_sex$gene <- gsub("\\..*", "", rownames(deg_sex))


####################################################################################

deg_by_clust = readRDS("Results/deg_by_clust_table_DecontX_final_DEG.All.RDS")

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

for (celltype in unique(deg_by_clust$cluster)){
  degs = deg_by_clust[deg_by_clust$cluster == celltype,]
  table(degs$change)
  x_range = max(abs(degs$avg_log2FC)) * 1.2
  
  ggplot(degs, aes(x = avg_log2FC, y = -log10(p_val_adj))) + 
    geom_point(aes(fill = change), 
               shape = 21, alpha = 0.75,
               na.rm = F, stroke = 0, size=2.5) + 
    # theme_classic(base_size = 15) + 
    theme_classic() +
    labs(title = paste("DEGs of", celltype ,"- Test vs Control", sep = " ")) + # Add a title
    xlab(expression(log[2]("Expr A" / "B"))) + # x-axis label
    ylab(expression(-log[10]("FDR"))) + # y-axis label
    geom_hline(yintercept = -log10(0.05), colour = "grey60", size = 1) + # Add FDR cutoffs
    #  geom_vline(xintercept = 0.25, linetype='dashed', colour = "grey40", size = 0.7) + # Add up-DEG cutoffs 
    geom_vline(xintercept = 0, colour = "grey80", size = 1.2) + # Add 0 lines
    #  geom_vline(xintercept = -0.25, linetype='dashed', colour = "grey40", size = 0.7) + # Add dn-DEG cutoffs 
    ggplot2::xlim(-x_range, x_range) + 
    geom_text_repel(d = subset(degs, change != "None"),
                    aes(label = gene),
                    size = 3.5, max.overlaps = 5) +
    theme(
      panel.grid.major = element_line(colour = "grey80", size = 0.5),
      panel.grid.minor = element_line(colour = "grey90", size = 0.25),
      # axis.text = element_text(size = 10),
      # axis.title = element_text(size = 12),
      # legend.text = element_text(size = 10),
      # legend.title = element_text(size =12),
      # title = element_text(size = 12),
      legend.background = element_rect(fill = "white"),
      legend.margin=margin(0.1,0.1,0.1,0.1)) +
    scale_fill_manual(values = c("Up" = "#E64B35", 
                                 "Down" = "#3182bd", 
                                 "None" = "grey"))
  ggsave(last_plot(), file = paste("Figures_Volcano/DEG_Volcano.", celltype, ".Test_vs_Control.pdf", sep = ""), width= 6, height = 5)
}



