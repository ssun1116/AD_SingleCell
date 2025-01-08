rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(dplyr)
library(SCENIC) 
library(AUCell)
library(RcisTarget)
library(SCopeLoomR) 
library(KernSmooth)
library(BiocParallel)
library(ggplot2)
library(data.table) ## may need re download.
library(grid)
library(ComplexHeatmap)
library(SingleCellExperiment)
library(Seurat)
set.seed(123)

### Species-specific database
#dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
#             "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")
# mc9nr: Motif collection version 9: 24k motifs
# for(featherURL in dbFiles){
#   download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
# }

## Load data
combined = readRDS("Data/data_combined_DecontX_Annotated_final.RDS")
combined = combined[!row.names(combined) %in% c("Malat1", "Gm26917", "Gm42418", "AY036118")] 


DimPlot(combined, label = T)
colnames(combined@meta.data)
table(combined$cond)

exprMat <- combined@assays$RNA@counts %>% as.matrix()
combined$celltype_cond = paste(combined$celltype, combined$cond.sex, sep = "_")
head(combined$celltype_cond)
cellInfo <- data.frame(seuratCluster = combined$celltype_cond)
loom <- build_loom("Combined_SCENIC_0513.loom", dgem = exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

# ### Initialize SCENIC settings
org <- "mgi" # or hgnc, or mgi
dbDir <- "Resources/cisTarget_databases/" # RcisTarget databases location
motifAnnotations_mgi = motifAnnotations
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, nCores=16) # Run on linux terminal.
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

loom <- open_loom("Combined_SCENIC_0513.loom")
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
#runGenie3(exprMat_filtered_log, scenicOptions)
exportsForArboreto(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
GRNBoost_out <- read.table("int/grn_output.tsv")
colnames(GRNBoost_out) <- c("TF", "Target", "weight")
saveRDS(GRNBoost_out,"int/1.4_GENIE3_linkList.Rds")

exprMat_log <- log2(exprMat+1)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat_log, skipHeatmap = TRUE)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

regulonAUC = readRDS("int/3.4_regulonAUC.Rds")
scenicOptions = readRDS("int/scenicOptions.Rds")
loom <- open_loom("Combined_SCENIC_0513.loom")
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster), function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

regulonActivity_highconf = regulonActivity_byCellType_Scaled[!grepl("extended", rownames(regulonActivity_byCellType_Scaled)),]

#ComplexHeatmap::Heatmap(t(regulonActivity_byCellType_Scaled), name="Regulon activity")
ComplexHeatmap::Heatmap(t(regulonActivity_highconf), name="Regulon activity")




### Exploring output
# Check files in folder 'output'
# Browse the output .loom file @ http://scope.aertslab.org

# tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
# aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
# 
# # Show TF expression:
# par(mfrow=c(2,3))
# AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat_log, 
#                         aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[7:12],],
#                         plots="Expression")

# # output/Step2_MotifEnrichment_preview.html in detail/subset:
# motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
# tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Sox8"]
# viewMotifs(tableSubset)
# 
# # output/Step2_regulonTargetsInfo.tsv in detail:
# regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
# tableSubset <- regulonTargetsInfo[TF=="Stat6" & highConfAnnot==TRUE]
# viewMotifs(tableSubset)
# 
# # Cell-type specific regulators (RSS):
# regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
# rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "CellType"], )
# rssPlot <- plotRSS(rss)
# plotly::ggplotly(rssPlot$plot)

seurat.cog = subset(seurat.all, subset = orig.ident == "chromium034")
exprMat <- seurat.cog@assays$RNA@counts %>% as.matrix()
cellInfo <- data.frame(seuratCluster = seurat.cog$clusterbysample)
loom <- build_loom("Seurat_Cog.Annoated.Final_1026.loom", dgem = exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

# ### Initialize SCENIC settings
org <- "mgi" # or hgnc, or mgi
dbDir <- "~/Dropbox/NGR_SNU_2019/scRNA_seq_Colab/Resources/cisTarget_databases/" # RcisTarget databases location
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, nCores=16) # Run on linux terminal.
saveRDS(scenicOptions, file="int/scenicOptions_Seurat.cog.Rds")

loom <- open_loom("Seurat_Cog.Annoated.Final_1026.loom")
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
#runGenie3(exprMat_filtered_log, scenicOptions)
exportsForArboreto(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
GRNBoost_out <- read.table("int/grn_output_Seurat.cog.tsv")
colnames(GRNBoost_out) <- c("TF", "Target", "weight")
saveRDS(GRNBoost_out,"int/1.4_GENIE3_linkList.Rds")

exprMat_log <- log2(exprMat+1)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat_log, skipHeatmap = TRUE)
saveRDS(scenicOptions, file="int/scenicOptions_Seurat.cog.Rds")

regulonAUC = readRDS("int/3.4_regulonAUC.Rds")

combined$celltype_cond = paste(combined$celltype, combined$cond, sep = "_")
head(combined$celltype_cond)
cellInfo <- data.frame(seuratCluster = combined$celltype_cond)

regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster), function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

regulonActivity_highconf = regulonActivity_byCellType_Scaled[!grepl("extended", rownames(regulonActivity_byCellType_Scaled)),]

ComplexHeatmap::Heatmap(t(regulonActivity_byCellType_Scaled), name="Regulon activity")
ComplexHeatmap::Heatmap(t(regulonActivity_highconf), name="Regulon activity")

