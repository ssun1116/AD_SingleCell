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

## Logcount mean - ANGPT2 grouping 
ad <- anndata::read_h5ad('Data/Pseudobulk_ROSMAP2_logcount_mean_0206.h5ad')
ad.raw = as.data.frame(ad$X)
ad.mtx = ad.raw %>% t()
ad.meta = as.data.frame(ad$obs)

ad.meta_endo = ad.meta[ad.meta$cluster_Allen2 == "Endothelial",]
ad.raw_endo = ad.raw[rownames(ad.meta_endo),]
ad.mtx_endo = ad.raw_endo %>% t()
ad.raw_endo_t = ad.mtx_endo %>% as.data.frame

## Quantile distribution for ANGPT2
quantile(ad.raw_endo$ANGPT2, probs = seq(0, 1, 0.25))

ad.raw_endo$ANGPT2_group = ifelse(ad.raw_endo$ANGPT2 > 1.3007325, "ANGPT2_High",
                                  ifelse(ad.raw_endo$ANGPT2 < 0.5207846, "ANGPT2_Low", "ANGPT2_Mid"))
table(ad.raw_endo$ANGPT2_group)
ad.raw_endo$ANGPT2_group = factor(ad.raw_endo$ANGPT2_group, 
                                  levels = rev(c("ANGPT2_Low", "ANGPT2_Mid", "ANGPT2_High")))
ad.raw_endo$sample = rownames(ad.raw_endo)
ad.raw_endo = ad.raw_endo[,c("ANGPT2_group", "sample")]
ad.raw_endo$sampleID = sapply(strsplit(as.character(ad.raw_endo$sample),'_'), "[", 1)
ad.raw_endo$sample <- NULL

## Raw count data
ad <- anndata::read_h5ad('Data/Pseudobulk_ROSMAP2_count_sum_0206.h5ad')
ad.raw = as.data.frame(ad$X)
ad.mtx = ad.raw %>% t()
ad.meta = as.data.frame(ad$obs)
ad.meta$rowname = rownames(ad.meta)

ad.meta_group = ad.meta[ad.meta$sampleID %in% ad.raw_endo$sampleID,]
ad.meta_group = merge(ad.raw_endo, ad.meta_group, by = "sampleID")
rownames(ad.meta_group) = ad.meta_group$rowname


## Comparison by ANGPT2 group
ad.meta_group$Age = as.numeric(levels(ad.meta_group$Age)[ad.meta_group$Age])
ad.meta_group$PMI = as.numeric(levels(ad.meta_group$PMI)[ad.meta_group$PMI])
ad.meta_group$age_death = as.numeric(levels(ad.meta_group$age_death)[ad.meta_group$age_death])
ad.meta_group$age_death = ifelse(is.na(ad.meta_group$age_death), 90, ad.meta_group$age_death)
ad.meta_group = ad.meta_group[!is.na(ad.meta_group$PMI),]

for (celltype in unique(ad.meta_group$cluster_Allen2)){
  if (celltype != "Others" & celltype != "SMC"){
    ad.meta_group_tmp = ad.meta_group[ad.meta_group$cluster_Allen2 == celltype & ad.meta_group$Diagnosis == "AD_dcfdx4",]
    
    ad.raw2 = ad.raw[rownames(ad.meta_group_tmp),]
    ad.raw2_t = ad.raw2 %>% t() %>% as.data.frame
    
    dds <- DESeqDataSetFromMatrix(ad.raw2_t,
                                  colData = ad.meta_group_tmp,
                                  design = ~ age_death + PMI + ANGPT2_group)
    
    dds$ANGPT2_group <- relevel(dds$ANGPT2_group, ref = "ANGPT2_Mid")
    
    dim(dds) # 35746 376
    # keep <- rowSums(counts(dds) >= 10) >= 4
    # table(keep)
    # dds <- dds[keep,]
    
    dds <- DESeq(dds)
    resultsNames(dds) # lists the coefficients
    res <- results(dds, name="ANGPT2_group_ANGPT2_High_vs_ANGPT2_Mid")
    # res <- lfcShrink(dds, coef="Diagnosis_AD_dcfdx4_vs_AD_dcfdx2", type="apeglm")
    
    res1 = data.frame(gene_name = rownames(res), baseMean=res$baseMean, log2FoldChange=res$log2FoldChange, 
                      lfcSE=res$lfcSE, pvalue=res$pvalue, padj=res$padj) 
    res1$Change = ifelse(res1$pvalue < 0.05 & abs(res1$log2FoldChange) > 0.25,
                         ifelse(res1$log2FoldChange > 0, "Up", "Down"), "None")
    res2 = res1[!is.na(res1$padj),] %>% arrange(-log2FoldChange)
    
    write_xlsx(res2, paste0("DESeq2_Pseudobulk_", celltype, "_AD_DV4_ANG2_High_vs_Mid_cov.age_death.pmi.xlsx"), 
               col_names = T)
  }
}



###########################################################
## Visualization
res.astro = read_excel("DESeq2_Pseudobulk_Endo_Diagnosis_cov.age_death.pmi.xlsx")

pp = fgsea::gmtPathways("GOBP_BLOOD_VESSEL_MORPHOGENESIS.v2023.2.Hs.gmt")
angio = read.table("GOBP_angiogenesis.txt", sep = "\t")


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


p


## CERAD

plot_title = "CERAD: Severe (Value 4) vs. Moderate (Value 3)"

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
  ylab(expression(-log[10]("FDR"))) + # y-axis label
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
