rm(list = ls())
library(dplyr)
library(tximport)
library(DESeq2)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(ggrastr)
library(Seurat)
library(writexl)
library(readxl)
library(anndata)
setwd("/data1/ROSMAP/")

## Analyze DEG (Raw count)
ad <- anndata::read_h5ad('Data/Pseudobulk_ROSMAP2_count_sum_0206.h5ad')
ad.raw = as.data.frame(ad$X)
ad.mtx = ad.raw %>% t()
ad.meta = as.data.frame(ad$obs)

ad.meta_endo = ad.meta[ad.meta$cluster_Allen2 == "Endothelial",]
ad.raw_endo = ad.raw[rownames(ad.meta_endo),]
ad.mtx_endo = ad.raw_endo %>% t()
ad.raw_endo_t = ad.mtx_endo %>% as.data.frame

## Comparison by Diagnosis (4 vs. 2)

ad.meta_endo.diag = ad.meta_endo[ad.meta_endo$Diagnosis %in% c("AD_dcfdx4", "AD_dcfdx2"),]
ad.raw_endo.diag = ad.raw_endo[rownames(ad.meta_endo.diag),]
ad.raw_endo.diag_t = ad.raw_endo.diag %>% t() %>% as.data.frame

ad.meta_endo.diag$Age = as.numeric(levels(ad.meta_endo.diag$Age)[ad.meta_endo.diag$Age])
ad.meta_endo.diag$PMI = as.numeric(levels(ad.meta_endo.diag$PMI)[ad.meta_endo.diag$PMI])
ad.meta_endo.diag$age_death = as.numeric(levels(ad.meta_endo.diag$age_death)[ad.meta_endo.diag$age_death])
ad.meta_endo.diag$age_death = ifelse(is.na(ad.meta_endo.diag$age_death), 90, ad.meta_endo.diag$age_death)

dds <- DESeqDataSetFromMatrix(ad.raw_endo.diag_t,
                                colData = ad.meta_endo.diag,
                                design = ~ age_death + PMI + Diagnosis)

dim(dds) # 35746 223
# keep <- rowSums(counts(dds) >= 10) >= 4
# table(keep)
# dds <- dds[keep,]

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="Diagnosis_AD_dcfdx4_vs_AD_dcfdx2")
# res <- lfcShrink(dds, coef="Diagnosis_AD_dcfdx4_vs_AD_dcfdx2", type="apeglm")

res1 = data.frame(gene_name = rownames(res), baseMean=res$baseMean, log2FoldChange=res$log2FoldChange, lfcSE=res$lfcSE, pvalue=res$pvalue, padj=res$padj) 
res1[res1$gene_name == "ANGPT2",]$log2FoldChange
res1[res1$gene_name == "ANGPT2",]$padj
res1[res1$gene_name == "ANGPT2",]$pvalue

res1$Change = ifelse(res1$pvalue < 0.05 & abs(res1$log2FoldChange) > 0.25,
                     ifelse(res1$log2FoldChange > 0, "Up", "Down"), "None")
res2 = res1[!is.na(res1$padj),] %>% arrange(-log2FoldChange)

write_xlsx(res2, "DESeq2_Pseudobulk_Endo_Diagnosis_cov.age_death.pmi.xlsx", col_names = T)



## Comparison by CERAD (Severe vs. Mild)

ad.meta_endo.cerad = ad.meta_endo.diag[ad.meta_endo.diag$Diagnosis == "AD_dcfdx4" & 
                                         ad.meta_endo.diag$ceradsc %in% c("3.0", "2.0", "1.0"),]

ad.meta_endo.cerad = ad.meta_endo.cerad[!is.na(ad.meta_endo.cerad$PMI), ] # 283

ad.raw_endo.cerad = ad.raw_endo[rownames(ad.meta_endo.cerad),]
ad.raw_endo.cerad_t = ad.raw_endo.cerad %>% t() %>% as.data.frame

#ad.meta_endo.cerad$Age = as.numeric(levels(ad.meta_endo.cerad$Age)[ad.meta_endo.cerad$Age])
#ad.meta_endo.cerad$PMI = as.numeric(levels(ad.meta_endo.cerad$PMI)[ad.meta_endo.cerad$PMI])
#ad.meta_endo.cerad$age_death = as.numeric(levels(ad.meta_endo.cerad$age_death)[ad.meta_endo.cerad$age_death])
#ad.meta_endo.cerad$age_death = ifelse(is.na(ad.meta_endo.cerad$age_death), 90, ad.meta_endo.cerad$age_death)

dds <- DESeqDataSetFromMatrix(ad.raw_endo.cerad_t,
                              colData = ad.meta_endo.cerad,
                              design = ~ age_death + PMI + ceradsc)
dds$ceradsc <- relevel(dds$ceradsc, ref = "3.0")

dim(dds) # 35746 116
# keep <- rowSums(counts(dds) >= 10) >= 4
# table(keep)
# dds <- dds[keep,]

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="ceradsc_2.0_vs_3.0")
# res <- lfcShrink(dds, coef="Diagnosis_AD_dcfdx4_vs_AD_dcfdx2", type="apeglm")

res1 = data.frame(gene_name = rownames(res), baseMean=res$baseMean, log2FoldChange=res$log2FoldChange, lfcSE=res$lfcSE, pvalue=res$pvalue, padj=res$padj) 
res1[res1$gene_name == "ANGPT2",]$log2FoldChange
res1[res1$gene_name == "ANGPT2",]$padj
res1[res1$gene_name == "ANGPT2",]$pvalue

res1$Change = ifelse(res1$pvalue < 0.05 & abs(res1$log2FoldChange) > 0.25,
                     ifelse(res1$log2FoldChange > 0, "Up", "Down"), "None")
res2 = res1[!is.na(res1$padj),] %>% arrange(-log2FoldChange)

write_xlsx(res2, "DESeq2_Pseudobulk_Endo_CERAD_cov.age_death.pmi.xlsx", col_names = T)


###########################################################
## Visualization
res.diag = read_excel("DESeq2_Pseudobulk_Endo_Diagnosis_cov.age_death.pmi.xlsx")
res.cog = read_excel("DESeq2_Pseudobulk_Endo_CERAD_cov.age_death.pmi.xlsx")

pp = fgsea::gmtPathways("GOBP_BLOOD_VESSEL_MORPHOGENESIS.v2023.2.Hs.gmt")
# angio = read.table("GOBP_angiogenesis.txt", sep = "\t")


## Diagnosis

plot_title = "Diagnosis: AD (Value 4) vs. MCI (Value 2)"

up_genes = res.diag %>% filter(Change != "None" & log2FoldChange > 0) %>% top_n(10, -log10(pvalue)) %>% pull(gene_name)
dn_genes = res.diag %>% filter(Change != "None" & log2FoldChange < 0) %>% top_n(10, -log10(pvalue)) %>% pull(gene_name)

up_genes = res.diag %>% filter(Change != "None" & log2FoldChange > 0) %>% pull(gene_name) %>% intersect(pp$GOBP_BLOOD_VESSEL_MORPHOGENESIS,.)
dn_genes = res.diag %>% filter(Change != "None" & log2FoldChange < 0) %>% pull(gene_name) %>% intersect(pp$GOBP_BLOOD_VESSEL_MORPHOGENESIS,.)

up_genes = res.diag %>% filter(Change != "None" & log2FoldChange > 0) %>% pull(gene_name) %>% intersect(angio$V1,.)
dn_genes = res.diag %>% filter(Change != "None" & log2FoldChange < 0) %>% pull(gene_name) %>% intersect(angio$V1,.)


#res.genes = res.diag %>% filter(Change != "None") %>% pull(gene_name)
#res.genes = c("ANGPT2", "ANGPT2", "TIE2", "TEK")

x_range = max(abs(res.diag$log2FoldChange)) * 1.5
p <- ggplot(res.diag, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point_rast(aes(fill = Change),
                  shape = 21, alpha = 0.75,
                  na.rm = F, stroke = 0, size=2.5) + # Make dots bigger
  theme_classic(base_size = 15) + # change theme
  labs(title = plot_title,) + # Add a title
  xlab(expression(log[2]("Expr A" / "B"))) + # x-axis label
  ylab(expression(-log[10]("FDR"))) + # y-axis label
  geom_hline(yintercept = 1.3, colour = "grey40", linetype='dashed', size = 0.7) + # Add cutoffs
  geom_vline(xintercept = 0.25, linetype='dashed', colour = "grey40", size = 0.7) + # Add 0 lines
  geom_vline(xintercept = -0.25, linetype='dashed', colour = "grey40", size = 0.7) + # Add 0 lines
  ggplot2::xlim(-x_range, x_range) +
  ggplot2::ylim(0, 1.1*max(-log10(res.diag$pvalue))) +
  # geom_text_repel(d = subset(res.diag, res.diag$gene_name %in% res.genes),
  #                 aes(label = gene_name),
  #                 size = 3.5, max.overlaps = 20) +
  geom_label_repel(
    data          = subset(res.diag, res.diag$gene_name %in% up_genes),
    nudge_x       = x_range*0.9 - subset(res.diag, res.diag$gene_name %in% up_genes)$log2FoldChange,
    segment.size  = 0.5,
    segment.color = "grey50",
    direction     = "y",
    aes(label = gene_name),
    size = 4, box.padding = 0.5, max.overlaps = Inf, fill = "white"
  ) +
  geom_label_repel(
    data          = subset(res.diag, res.diag$gene_name %in% dn_genes),
    nudge_x       = -x_range*0.9 - subset(res.diag, res.diag$gene_name %in% dn_genes)$log2FoldChange,
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


ggsave(p, "")


## CERAD

plot_title = "CERAD: Severe (Value 1) vs. Mild (Value 3)"

up_genes = res.cog %>% filter(Change != "None" & log2FoldChange > 0) %>% top_n(10, -log10(pvalue)) %>% pull(gene_name)
dn_genes = res.cog %>% filter(Change != "None" & log2FoldChange < 0) %>% top_n(10, -log10(pvalue)) %>% pull(gene_name)

up_genes = res.cog %>% filter(Change != "None" & log2FoldChange > 0) %>% pull(gene_name) %>% intersect(pp$GOBP_BLOOD_VESSEL_MORPHOGENESIS,.)
dn_genes = res.cog %>% filter(Change != "None" & log2FoldChange < 0) %>% pull(gene_name) %>% intersect(pp$GOBP_BLOOD_VESSEL_MORPHOGENESIS,.)


#res.genes = res.cog %>% filter(Change != "None") %>% pull(gene_name)
#res.genes = c("ANGPT2", "ANGPT2", "TIE2", "TEK")

x_range = max(abs(res.cog$log2FoldChange)) * 1.5
p <- ggplot(res.cog, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point_rast(aes(fill = Change),
                  shape = 21, alpha = 0.75,
                  na.rm = F, stroke = 0, size=2.5) + # Make dots bigger
  theme_classic(base_size = 15) + # change theme
  labs(title = plot_title,) + # Add a title
  xlab(expression(log[2]("Expr A" / "B"))) + # x-axis label
  ylab(expression(-log[10]("p-value"))) + # y-axis label
  geom_hline(yintercept = 1.3, colour = "grey40", linetype='dashed', size = 0.7) + # Add cutoffs
  geom_vline(xintercept = 0.25, linetype='dashed', colour = "grey40", size = 0.7) + # Add 0 lines
  geom_vline(xintercept = -0.25, linetype='dashed', colour = "grey40", size = 0.7) + # Add 0 lines
  ggplot2::xlim(-x_range, x_range) +
  ggplot2::ylim(0, 1.1*max(-log10(res.cog$pvalue))) +
  # geom_text_repel(d = subset(res.cog, res.cog$gene_name %in% res.genes),
  #                 aes(label = gene_name),
  #                 size = 3.5, max.overlaps = 20) +
  geom_label_repel(
    data          = subset(res.cog, res.cog$gene_name %in% up_genes),
    nudge_x       = x_range*0.9 - subset(res.cog, res.cog$gene_name %in% up_genes)$log2FoldChange,
    segment.size  = 0.5,
    segment.color = "grey50",
    direction     = "y",
    aes(label = gene_name),
    size = 4, box.padding = 0.5, max.overlaps = Inf, fill = "white"
  ) +
  geom_label_repel(
    data          = subset(res.cog, res.cog$gene_name %in% dn_genes),
    nudge_x       = -x_range*0.9 - subset(res.cog, res.cog$gene_name %in% dn_genes)$log2FoldChange,
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


p
