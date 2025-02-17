rm(list = ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(writexl)

combined = readRDS("Data/Seurat_Integrated_DecontX_Edit_DraftAnnot_0219.RDS")

# Multiplot code
plot_multi <- function(obj = obj, gene = gene, group = "final_annot"){
  p1 <- FeaturePlot(obj, features = gene, label = T)
  # obj2 = subset(obj, subset = celltype != "Low-quality")
  p2 <- VlnPlot(obj, features = gene, pt.size = 0) +
    theme(axis.title.x = element_blank(),
          legend.position = "None",
          plot.title = element_blank())
  p = plot_grid(p1, p2, ncol = 1, rel_heights = c(1.75,1))
  
  ggsave(paste0("Figures/GeneExpr_", gene, "_Multiplot.pdf"), p, width = 4.5, height = 6.5)
}

# Jeungeun-curated marker
runtype = "EL"

geneset = c("Slc17a7", "Camk2a", "Slc17a6", "Wfs1", "Map3k15", "Crlf1", "Nts", "Gad1", "Pvalb", "Lamp5")
geneset2 = c("Gad2", "Slc32a1")
geneset3 = c("Plekhg1", "Lefty1", "Gpc3", "Dlk1", "Ramp3", "Map3k15", 
             "Iyd", "Prox1", "Trp73")
geneset4 = c("Angpt2", "Pecam1", "Cdh5", "Angptl1", "Tek")

for (gene in geneset4){
  plot_multi(combined, gene, "final_annot")
}






