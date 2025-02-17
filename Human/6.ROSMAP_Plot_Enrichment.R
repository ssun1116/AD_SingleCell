rm(list = ls())
library(ggplot2)
library(readxl)
library(cowplot)
library(ggplot2)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x 
}

df_lin1 = read_excel("enrichment_result.xlsx", sheet = 2)
df_lin2 = read_excel("enrichment_result.xlsx", sheet = 7)


df_lin1 = read_excel("Enrichment_DESeq2_Pseudobulk_Endothelial_Diagnosis_cov.age_death.pmi.xlsx", sheet = 1)
df_lin2 = read_excel("Enrichment_DESeq2_Pseudobulk_Endothelial_Diagnosis_cov.age_death.pmi.xlsx", sheet = 2)


lin1 = c(
  "regulation of membrane potential",
  "cell junction organization",
  "neuron projection morphogenesis",
  "enzyme-linked receptor protein signaling pathway",
  "response to growth factor",
  "regulation of transmembrane transport ",
  "response to hypoxia"
)
lin2 = c(
  "response to cytokine",
  "response to virus",
  "wound healing",
  "cytokine-mediated signaling pathway",
  "positive regulation of cellular component movement",
  "cellular oxidant detoxification",
  "macrophage activation"
)

## lineage
df_lin1_select = df_lin1[df_lin1$Description %in% lin1,]
df_lin1_select$Description_ID = paste0(firstup(df_lin1_select$Description), "\n(",  df_lin1_select$ID, ")")

p1 = ggplot(df_lin1_select, aes(x = reorder(Description_ID, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", 
           color = 'black', lwd = 0.3,
           fill = '#C593C2', 
           alpha = 1, position = "identity", show.legend = F) +
  labs(x = '', y = '-log10(FDR)') +
  coord_flip() +
  #  coord_fixed(ratio = 0.5) +
  theme_bw(base_size = 7) +
  # ylim(0, 40) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x=element_text(size = 5, colour="black", angle = 0, vjust = 3, hjust = 0),
    axis.text.y = element_text(size = 5, colour = "black"),
    axis.title.x = element_text(size = 5, vjust = 0.1),
    plot.title = element_text(size = 6, face = "bold", vjust = 1.2, hjust = 0.5),
    legend.position="none",
    axis.ticks = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"))



## lineage 1 late

df_lin2_select = df_lin2[df_lin2$Description %in% lin2,]
df_lin2_select$Description_ID = paste0(firstup(df_lin2_select$Description), "\n(",  df_lin2_select$ID, ")")

p2 = ggplot(df_lin2_select, aes(x = reorder(Description_ID, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", 
           color = 'black', lwd = 0.3,
           fill = '#90C0DF', 
           alpha = 1, position = "identity", show.legend = F) +
  labs(x = '', y = '-log10(FDR)') +
  coord_flip() +
  #  coord_fixed(ratio = 0.5) +
  theme_bw(base_size = 7) +
  # ylim(0, 40) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x=element_text(size = 5, colour="black", angle = 0, vjust = 3, hjust = 0),
    axis.text.y = element_text(size = 5, colour = "black"),
    axis.title.x = element_text(size = 5, vjust = 0.1),
    plot.title = element_text(size = 6, face = "bold", vjust = 1.2, hjust = 0.5),
    legend.position="none",
    axis.ticks = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"))




## lineage 2 early

df_lin3_select = df_lin3[df_lin3$Description %in% lin3,]
df_lin3_select$Description_ID = paste0(firstup(df_lin3_select$Description), "\n(",  df_lin3_select$ID, ")")

p3 = ggplot(df_lin3_select, aes(x = reorder(Description_ID, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", 
           color = 'black', lwd = 0.3,
           fill = '#C593C2', 
           alpha = 1, position = "identity", show.legend = F) +
  labs(x = '', y = '-log10(FDR)') +
  coord_flip() +
  #  coord_fixed(ratio = 0.5) +
  theme_bw(base_size = 7) +
  # ylim(0, 40) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x=element_text(size = 5, colour="black", angle = 0, vjust = 3, hjust = 0),
    axis.text.y = element_text(size = 5, colour = "black"),
    axis.title.x = element_text(size = 5, vjust = 0.1),
    plot.title = element_text(size = 6, face = "bold", vjust = 1.2, hjust = 0.5),
    legend.position="none",
    axis.ticks = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"))



p = plot_grid(p1, p2,  nrow = 1)
ggsave(p, filename = "GO_Analysis.pdf", width = 10, height = 3)

########

df_lin1 = read_excel("GO_Analysis_Lineages_All_under_20_lineage2_time_wilcoxon_1000.xlsx", sheet = 1)
df_lin2 = read_excel("GO_Analysis_Lineages_All_under_20_lineage2_time_wilcoxon_1000.xlsx", sheet = 2)
df_lin3 = read_excel("/GO_Analysis_Lineages_All_under_20_lineage2_time_wilcoxon_1000.xlsx", sheet = 3)

lin1 = c(
  "chemical synaptic transmission",
  "synaptic signaling",
  "neuron projection morphogenesis",
  "neuron projection development",
  "synaptic transmission, glutamatergic",
  "synapse organization",
  "neurotransmitter secretion",
  "cytoskeletal protein binding",
  "regulation of cytoskeleton organization",
  "regulation of actin cytoskeleton organization"
)
lin2 = c(
  "chromosome organization",
  "mitotic cell cycle process",
  "regulation of cell cycle process",
  "chromosomal region",
  "translation",
  "chromatin organization",
  "chromatin remodeling",
  "regulation of canonical Wnt signaling pathway",
  "tissue morphogenesis",
  "regulation of Wnt signaling pathway"
)
lin3 = c(
  "chromosome organization",
  "spindle assembly",
  "translation",
  "chromosome segregation",
  "mitotic cell cycle",
  "microtubule cytoskeleton organization",
  "microtubule-based process",
  "nuclear chromosome",
  "microtubule",
  "regulation of cell cycle process"
)


## lineage
df_lin1_select = df_lin1[df_lin1$Description %in% lin1,]
df_lin1_select$Description_ID = paste0(firstup(df_lin1_select$Description), "\n(",  df_lin1_select$ID, ")")

p1 = ggplot(df_lin1_select, aes(x = reorder(Description_ID, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", 
           color = 'black', lwd = 0.3,
           fill = '#ff7f0eff', 
           alpha = 1, position = "identity", show.legend = F) +
  labs(x = '', y = '-log10(FDR)') +
  coord_flip() +
  #  coord_fixed(ratio = 0.5) +
  theme_bw(base_size = 7) +
  # ylim(0, 40) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x=element_text(size = 5, colour="black", angle = 0, vjust = 3, hjust = 0),
    axis.text.y = element_text(size = 5, colour = "black"),
    axis.title.x = element_text(size = 5, vjust = 0.1),
    plot.title = element_text(size = 6, face = "bold", vjust = 1.2, hjust = 0.5),
    legend.position="none",
    axis.ticks = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"))



## lineage 1 late

df_lin2_select = df_lin2[df_lin2$Description %in% lin2,]
df_lin2_select$Description_ID = paste0(firstup(df_lin2_select$Description), "\n(",  df_lin2_select$ID, ")")

p2 = ggplot(df_lin2_select, aes(x = reorder(Description_ID, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", 
           color = 'black', lwd = 0.3,
           fill = '#ff7f0eff', 
           alpha = 0.7, position = "identity", show.legend = F) +
  labs(x = '', y = '-log10(FDR)') +
  coord_flip() +
  #  coord_fixed(ratio = 0.5) +
  theme_bw(base_size = 7) +
  # ylim(0, 40) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x=element_text(size = 5, colour="black", angle = 0, vjust = 3, hjust = 0),
    axis.text.y = element_text(size = 5, colour = "black"),
    axis.title.x = element_text(size = 5, vjust = 0.1),
    plot.title = element_text(size = 6, face = "bold", vjust = 1.2, hjust = 0.5),
    legend.position="none",
    axis.ticks = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"))




## lineage 2 early

df_lin3_select = df_lin3[df_lin3$Description %in% lin3,]
df_lin3_select$Description_ID = paste0(firstup(df_lin3_select$Description), "\n(",  df_lin3_select$ID, ")")

p3 = ggplot(df_lin3_select, aes(x = reorder(Description_ID, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", 
           color = 'black', lwd = 0.3,
           fill = '#ff7f0eff', 
           alpha = 0.2, position = "identity", show.legend = F) +
  labs(x = '', y = '-log10(FDR)') +
  coord_flip() +
  #  coord_fixed(ratio = 0.5) +
  theme_bw(base_size = 7) +
  # ylim(0, 40) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    axis.text.x=element_text(size = 5, colour="black", angle = 0, vjust = 3, hjust = 0),
    axis.text.y = element_text(size = 5, colour = "black"),
    axis.title.x = element_text(size = 5, vjust = 0.1),
    plot.title = element_text(size = 6, face = "bold", vjust = 1.2, hjust = 0.5),
    legend.position="none",
    axis.ticks = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"))



p = plot_grid(p3, p2, p1, nrow = 1)
ggsave(p, filename = "GO_Analysis_Lineages_All_under_20_lineage2_time_wilcoxon_1000.pdf", width = 10, height = 2)



