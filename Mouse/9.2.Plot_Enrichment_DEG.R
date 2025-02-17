rm(list = ls())
options(stringsAsFactors = F)
library(ggplot2)
library(dplyr)
library(cowplot)
library(readxl)
library(writexl)

## Tables
d1_up = read_excel("Tables/Share/cluster_Ext (Cortex)_enrichment_result_EL.xlsx", sheet = 1)
d1_dn = read_excel("Tables/Share/cluster_Ext (Cortex)_enrichment_result_EL.xlsx", sheet = 2)

d2_up = read_excel("Tables/Share/cluster_Astro_enrichment_result_EL.xlsx", sheet = 1)
d2_dn = read_excel("Tables/Share/cluster_Astro_enrichment_result_EL.xlsx", sheet = 2)

d3_up = read_excel("Tables/Share/cluster_Micro_enrichment_result_EL.xlsx", sheet = 1)
d3_dn = read_excel("Tables/Share/cluster_Micro_enrichment_result_EL.xlsx", sheet = 2)

## Terms
d1_up_term = c("neuron projection morphogenesis", "cell morphogenesis involved in neuron differentiation", "synaptic signaling", 
               "synapse organization", "dendrite morphogenesis", "modulation of chemical synaptic transmission", 
               "regulation of trans-synaptic signaling", "synaptic transmission, glutamatergic", "regulation of ion transport",
               "negative regulation of cell migration", "negative regulation of cellular component movement")
  
d1_dn_term = c("synaptic signaling", "chemical synaptic transmission", "anterograde trans-synaptic signaling", "trans-synaptic signaling", 
               "behavior", "learning or memory", "modulation of chemical synaptic transmission", "regulation of trans-synaptic signaling", 
               "regulation of neuron projection development", "cognition", "learning", 
               "synapse organization", "memory", "positive regulation of nervous system development", "regulation of axonogenesis")  

d2_up_term = c("synapse organization", "synapse assembly", "cell-cell adhesion via plasma-membrane adhesion molecules", 
               "cell junction assembly", "synaptic signaling", "central nervous system neuron differentiation", "neuron migration")

d2_dn_term = c("negative regulation of apoptotic process", "glial cell development", "negative regulation of neuron apoptotic process", 
               "glial cell differentiation", "negative regulation of neuron death", "gliogenesis", "MAPK cascade")  

d3_up_term = c("actin filament-based process", "actin cytoskeleton organization", "regulation of cell projection organization", 
               "cell junction organization", "cell morphogenesis involved in differentiation", "neuron projection morphogenesis", 
               "plasma membrane bounded cell projection morphogenesis", "cell projection morphogenesis")

d3_dn_term = c("positive regulation of immune response", "positive regulation of response to biotic stimulus", "innate immune response", 
               "regulation of adaptive immune response", "positive regulation of response to external stimulus", 
               "positive regulation of innate immune response", "cellular response to molecule of bacterial origin", "response to lipopolysaccharide")  


## Subset tables
d1_up_sub = d1_up[d1_up$Description %in% d1_up_term,] # 11
d1_dn_sub = d1_dn[d1_dn$Description %in% d1_dn_term,] # 15

d2_up_sub = d2_up[d2_up$Description %in% d2_up_term,] # 7
d2_dn_sub = d2_dn[d2_dn$Description %in% d2_dn_term,] # 7

d3_up_sub = d3_up[d3_up$Description %in% d3_up_term,] # 8
d3_dn_sub = d3_dn[d3_dn$Description %in% d3_dn_term,] # 8


## Plot
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x 
}

## Ext cortex
d1_up_sub$Description_ID = paste0(firstup(d1_up_sub$Description), "\n(",  d1_up_sub$ID, ")")

p1.1 = ggplot(d1_up_sub, aes(x = reorder(Description_ID, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", 
           color = 'black', lwd = 0.3,
           fill = '#C593C2', 
           alpha = 1, position = "identity", show.legend = F) +
  labs(x = '', y = '-log10(FDR)', title = "Up-regulated DEGs in Ext (Cortex)") +
  coord_flip() +
  #  coord_fixed(ratio = 0.5) +
  theme_bw(base_size = 10) +
  # ylim(0, 40) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x=element_text(size = 10, colour="black", angle = 0, vjust = 2, hjust = 0),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 10, vjust = 0.1),
    plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
    legend.position="none",
    axis.ticks = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

d1_dn_sub$Description_ID = paste0(firstup(d1_dn_sub$Description), "\n(",  d1_dn_sub$ID, ")")

p1.2 = ggplot(d1_dn_sub, aes(x = reorder(Description_ID, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", 
           color = 'black', lwd = 0.3,
           fill = '#90C0DF', 
           alpha = 1, position = "identity", show.legend = F) +
  labs(x = '', y = '-log10(FDR)', title = "Down-regulated DEGs in Ext (Cortex)") +
  coord_flip() +
  #  coord_fixed(ratio = 0.5) +
  theme_bw(base_size = 10) +
  # ylim(0, 40) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x=element_text(size = 10, colour="black", angle = 0, vjust = 2, hjust = 0),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 10, vjust = 0.1),
    plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
    legend.position="none",
    axis.ticks = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

p1 = plot_grid(p1.1, p1.2, nrow = 1)
ggsave(plot = p1, filename = "Figures/Plot_cluster_Ext (Cortex)_enrichment_result_EL.pdf", width = 15, height = 5)

## Astrocyte

d2_up_sub$Description_ID = paste0(firstup(d2_up_sub$Description), "\n(",  d2_up_sub$ID, ")")

p2.1 = ggplot(d2_up_sub, aes(x = reorder(Description_ID, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", 
           color = 'black', lwd = 0.3,
           fill = '#C593C2', 
           alpha = 1, position = "identity", show.legend = F) +
  labs(x = '', y = '-log10(FDR)', title = "Up-regulated DEGs in Astrocyte") +
  coord_flip() +
  #  coord_fixed(ratio = 0.5) +
  theme_bw(base_size = 10) +
  # ylim(0, 40) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x=element_text(size = 10, colour="black", angle = 0, vjust = 2, hjust = 0),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 10, vjust = 0.1),
    plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
    legend.position="none",
    axis.ticks = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

d2_dn_sub$Description_ID = paste0(firstup(d2_dn_sub$Description), "\n(",  d2_dn_sub$ID, ")")

p2.2 = ggplot(d2_dn_sub, aes(x = reorder(Description_ID, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", 
           color = 'black', lwd = 0.3,
           fill = '#90C0DF', 
           alpha = 1, position = "identity", show.legend = F) +
  labs(x = '', y = '-log10(FDR)', title = "Down-regulated DEGs in Astrocyte") +
  coord_flip() +
  #  coord_fixed(ratio = 0.5) +
  theme_bw(base_size = 10) +
  # ylim(0, 40) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x=element_text(size = 10, colour="black", angle = 0, vjust = 2, hjust = 0),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 10, vjust = 0.1),
    plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
    legend.position="none",
    axis.ticks = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

p2 = plot_grid(p2.1, p2.2, nrow = 1)
ggsave(plot = p2, filename = "Figures/Plot_cluster_Astrocyte_enrichment_result_EL.pdf", width = 15, height = 4)

## Microglia

d3_up_sub$Description_ID = paste0(firstup(d3_up_sub$Description), "\n(",  d3_up_sub$ID, ")")

p3.1 = ggplot(d3_up_sub, aes(x = reorder(Description_ID, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", 
           color = 'black', lwd = 0.3,
           fill = '#C593C2', 
           alpha = 1, position = "identity", show.legend = F) +
  labs(x = '', y = '-log10(FDR)', title = "Up-regulated DEGs in Microglia") +
  coord_flip() +
  #  coord_fixed(ratio = 0.5) +
  theme_bw(base_size = 10) +
  # ylim(0, 40) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x=element_text(size = 10, colour="black", angle = 0, vjust = 2, hjust = 0),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 10, vjust = 0.1),
    plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
    legend.position="none",
    axis.ticks = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

d3_dn_sub$Description_ID = paste0(firstup(d3_dn_sub$Description), "\n(",  d3_dn_sub$ID, ")")

p3.2 = ggplot(d3_dn_sub, aes(x = reorder(Description_ID, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", 
           color = 'black', lwd = 0.3,
           fill = '#90C0DF', 
           alpha = 1, position = "identity", show.legend = F) +
  labs(x = '', y = '-log10(FDR)', title = "Down-regulated DEGs in Microglia") +
  coord_flip() +
  #  coord_fixed(ratio = 0.5) +
  theme_bw(base_size = 10) +
  # ylim(0, 40) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x=element_text(size = 10, colour="black", angle = 0, vjust = 2, hjust = 0),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 10, vjust = 0.1),
    plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
    legend.position="none",
    axis.ticks = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

p3 = plot_grid(p3.1, p3.2, nrow = 1)
ggsave(plot = p3, filename = "Figures/Plot_cluster_Microglia_enrichment_result_EL.pdf", width = 15, height = 4)



## 240812 Update.

## Tables
d1_up = read_excel("Tables/Enrichment_DESeq2_Pseudobulk_Endothelial_Diagnosis_cov.age_death.pmi EL MK final.xlsx", sheet = 1)
d1_dn = read_excel("Tables/Enrichment_DESeq2_Pseudobulk_Endothelial_Diagnosis_cov.age_death.pmi EL MK final.xlsx", sheet = 2)
d1_gse = read_excel("Tables/Enrichment_DESeq2_Pseudobulk_Endothelial_Diagnosis_cov.age_death.pmi EL MK final.xlsx", sheet = 3)

## Terms
d1_up_term = c("synaptic signaling", "modulation of chemical synaptic transmission", "cell junction organization", 
               "synapse organization", "cell junction assembly", "synapse assembly", "neuron projection morphogenesis")

d1_dn_term = c("response to cytokine", "response to virus", "cytokine-mediated signaling pathway", 
               "regulation of body fluid levels", "cellular oxidant detoxification", "lipid biosynthetic process",
               "macrophage activation", "cellular response to cytokine stimulus")  


## Subset tables
d1_up_sub = d1_up[d1_up$Description %in% d1_up_term,] # 7
d1_dn_sub = d1_dn[d1_dn$Description %in% d1_dn_term,] # 8


## Plot
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x 
}


d1_up_sub$Description_ID = paste0(firstup(d1_up_sub$Description), "\n(",  d1_up_sub$ID, ")")

p1.1 = ggplot(d1_up_sub, aes(x = reorder(Description_ID, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", 
           color = 'black', lwd = 0.3,
           fill = '#C593C2', 
           alpha = 1, position = "identity", show.legend = F) +
  labs(x = '', y = '-log10(FDR)', title = "Up-regulated DEGs") +
  coord_flip() +
  #  coord_fixed(ratio = 0.5) +
  theme_bw(base_size = 10) +
  # ylim(0, 40) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x=element_text(size = 10, colour="black", angle = 0, vjust = 2, hjust = 0),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 10, vjust = 0.1),
    plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
    legend.position="none",
    axis.ticks = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

d1_dn_sub$Description_ID = paste0(firstup(d1_dn_sub$Description), "\n(",  d1_dn_sub$ID, ")")

p1.2 = ggplot(d1_dn_sub, aes(x = reorder(Description_ID, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", 
           color = 'black', lwd = 0.3,
           fill = '#90C0DF', 
           alpha = 1, position = "identity", show.legend = F) +
  labs(x = '', y = '-log10(FDR)', title = "Down-regulated DEGs") +
  coord_flip() +
  #  coord_fixed(ratio = 0.5) +
  theme_bw(base_size = 10) +
  # ylim(0, 40) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x=element_text(size = 10, colour="black", angle = 0, vjust = 2, hjust = 0),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 10, vjust = 0.1),
    plot.title = element_text(size = 12, face = "bold", vjust = 1.2, hjust = 0.5),
    legend.position="none",
    axis.ticks = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

p1 = plot_grid(p1.1, p1.2, nrow = 1)
ggsave(plot = p1, filename = "Figures/Plot_Enrichment_DESeq2_Pseudobulk_Endothelial_Diagnosis_cov.age_death.pdf", width = 10, height = 3.5)


