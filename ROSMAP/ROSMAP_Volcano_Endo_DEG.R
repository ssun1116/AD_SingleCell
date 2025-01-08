rm(list = ls())
library(dplyr)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(ggrastr)
library(readxl)

res.diag = read_excel("Tables/DESeq2_Pseudobulk_Endo_Diagnosis_cov.age_death.pmi.xlsx")
res.cerad = read_excel("Tables/DESeq2_Pseudobulk_Endo_CERAD_cov.age_death.pmi.xlsx")

res2 = res.cerad

up_genes = res2 %>% filter(Change != "None" & log2FoldChange > 0) %>% top_n(10, -log10(pvalue)) %>% pull(gene_name)
dn_genes = res2 %>% filter(Change != "None" & log2FoldChange < 0) %>% top_n(10, -log10(pvalue)) %>% pull(gene_name)
#res.genes = res2 %>% filter(Change != "None") %>% pull(gene_name)
#res.genes = c("ANGPT2", "ANGPT2", "TIE2", "TEK")

perm.genes = c("TJP1", "PLVAP", "PDE2A", "GPR4", "YES1", "SRC", "PTP4A3", "OCLN", "CTNNBIP1", "BMP6", "FERMT2", "ZEB2", "C2CD4A", "VEGFA", "TJP2", "MYLK3", "ADORA2A", "NPR1", "TRPV4", "SLIT2", "CLDN5", "SH3GL2", "ANGPT1", "PLEC", "FGFBP3", "TGFB1", "C2CD4B", "AKAP12", "TJP3", "DDAH1", "NPPB", "BCR", "CDH5", "TACR2", "TACR1", "HRH1", "ADM", "TEK", "RAMP2", "ABCC8", "ARHGAP35", "PTPRJ", "APOE", "UCN", "CEACAM1", "CXCR2", "BDKRB2", "PDE3A", "AZU1", "ANGPT2")

vegf.genes = c("FLT3", "PRKD2", "PDGFRB", "VEGFD", "FLT1", "NRP2", "NUS1", "PIK3CA", "FBXW7-AS1", "VEGFA", "PGF", "VEGFC", "VEGFB", "MYO1C", "PIK3CD", "PDGFRA", "PRKD1", "HSPB1", "GAB1", "FOXC1", "NRP1", "RELA", "KDR", "FLT4", "MCAM", "NTN1" , "UNC5B", "ROBO4" )

vegf.reg.genes = c("SMOC2", "TNXB", "HRG", "DAB2IP", "SPRY2", "IL12B", "IL12A", "PTP4A3", "DCN", "PIK3CB", "ADAMTS3", "XDH", "CADM4", "ITGB1", "ROBO1", "VEGFA", "SEMA6A", "MYO1C", "CCBE1", "ADGRA2", "DLL1", "JCAD", "CD63", "EMILIN1" )

tie.genes = c("ANGPT2", "GRB2", "PIK3R2", "GRB14", "PIK3R1", "PIK3CB", "SOS1", "PTPN11", "HRAS", "PIK3CA", "SHC1", "KRAS", "ANGPT1", "GRB7", "DOK2", "NRAS", "TEK", "ANGPT4")

plot_title = "Diagnosis: AD (Value 4) vs. MCI (Value 2)"
plot_title = "Diagnosis: Regulation of Vascular Permeability"

plot_title = "CERAD: VEGF signaling pathway"
plot_title = "Diagnosis: Regulation of VEGF signaling pathway"
plot_title = "CERAD: TIE signaling pathway"

plot_title = "CERAD: Severe (Value 1.0) vs. Mild (Value 3.0)"
plot_title = "CERAD: Regulation of Vascular Permeability"

genes = perm.genes
x_range = max(abs(res2$log2FoldChange)) * 1.5
p <- ggplot(res2, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point_rast(aes(fill = Change),
                  shape = 21, alpha = 0.75,
                  na.rm = F, stroke = 0, size=2.5) + # Make dots bigger
  theme_classic(base_size = 15) + # change theme
  labs(title = plot_title,) + # Add a title
  xlab(expression(log[2]("Expr A" / "B"))) + # x-axis label
  ylab(expression(-log[10]("P-value"))) + # y-axis label
  geom_hline(yintercept = 1.3, colour = "grey40", linetype='dashed', size = 0.7) + # Add cutoffs
  geom_vline(xintercept = 0.25, linetype='dashed', colour = "grey40", size = 0.7) + # Add 0 lines
  geom_vline(xintercept = -0.25, linetype='dashed', colour = "grey40", size = 0.7) + # Add 0 lines
  ggplot2::xlim(-x_range, x_range) +
  ggplot2::ylim(0, 1.1*max(-log10(res2$pvalue))) +
  # geom_text_repel(d = subset(res2, res2$gene_name %in% perm.genes),
  #                 aes(label = gene_name),
  #                 size = 3.5, max.overlaps = 20) +
  geom_label_repel(
    data          = subset(res2, res2$Change == "Up" & res2$gene_name %in% genes),
    nudge_x       = x_range*0.9 - subset(res2, res2$Change == "Up" & res2$gene_name %in% genes)$log2FoldChange,
    segment.size  = 0.5,
    segment.color = "grey50",
    direction     = "y",
    aes(label = gene_name),
    size = 4, box.padding = 0.5, max.overlaps = Inf, fill = "white"
  ) +
  geom_label_repel(
    data          = subset(res2, res2$Change == "Down" & res2$gene_name %in% genes),
    nudge_x       = -x_range*0.9 - subset(res2, res2$Change == "Down" & res2$gene_name %in% genes)$log2FoldChange,
    segment.size  = 0.5,
    segment.color = "grey50",
    direction     = "y",
    aes(label = gene_name),
    size = 4, box.padding = 0.5, max.overlaps = Inf, fill = "white"
  ) +
theme(
  panel.grid.major = element_line(colour = "grey80", size = 0.5),
  panel.grid.minor = element_line(colour = "grey90", size = 0.25),
  axis.text = element_text(size = 10),
  axis.title = element_text(size = 12),
  legend.text = element_text(size = 10),
  legend.title = element_text(size =12),
  title = element_text(size = 12),
  # legend.position = c(0.5, 0.9),
  # legend.direction = "horizontal",
  legend.background = element_rect(fill = "white"),
  legend.margin=margin(0.1,0.1,0.1,0.1)) +
  scale_fill_manual(values = c("Up" = "#E64B35",
                               "Down" = "#3182bd",
                               "None" = "grey"))


ggsave(p, file = "Figures/Volcano_ROSMAP_CERAD_1_vs_3_Perm.pdf", 
       width = 6, height = 5)
