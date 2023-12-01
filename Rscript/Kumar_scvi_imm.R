library(tidyverse)
library(data.table)
library(Matrix)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratWrappers)
library(msigdbr)
library(fgsea)
library(speckle)
library(ggsignif)

# clean immune subset----
seu_imm <- readRDS("RDSfiles/seu_imm_scvi.RDS")

seu_imm <-DietSeurat(seu_imm) 
seu_imm <- JoinLayers(seu_imm)
seu_imm[["RNA"]]$data <- NULL
seu_imm[["RNA"]]$scale.data <- NULL
seu_imm[["RNA"]] <- split(seu_imm[["RNA"]], f = seu_imm$sample_ID)
seu_imm

seu_imm <- NormalizeData(seu_imm, verbose = FALSE)
seu_imm <- FindVariableFeatures(seu_imm, verbose = FALSE)
seu_imm <- ScaleData(seu_imm, verbose = FALSE)
seu_imm <- RunPCA(seu_imm, npcs = 30, verbose = FALSE)

seu_imm <- IntegrateLayers(
  object = seu_imm, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/Users/kinoshitahiroto/miniconda3/envs/scvi-env",
  verbose = TRUE
)

seu_imm <- FindNeighbors(seu_imm, reduction = "integrated.scvi", dims = 1:30, verbose = FALSE)
seu_imm <- FindClusters(seu_imm, resolution = 2, cluster.name = "scvi_clusters", verbose = FALSE)
seu_imm <- RunUMAP(seu_imm, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", verbose = FALSE)

DimPlot(seu_imm, label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu_imm, label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
DimPlot(seu_imm, group.by = "sample_ID") + NoAxes()

FeaturePlot(seu_imm, reduction = "umap.scvi", features = "EPCAM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "MUC5AC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "TFF1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "MUC6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "LIPF", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "PGA3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
# FeaturePlot(seu_imm, reduction = "umap.scvi", features = "REG4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
# FeaturePlot(seu_imm, reduction = "umap.scvi", features = "TFF3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
# FeaturePlot(seu_imm, reduction = "umap.scvi", features = "CHGA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
# FeaturePlot(seu_imm, reduction = "umap.scvi", features = "ATP4B", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "KRT18", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
# FeaturePlot(seu_imm, reduction = "umap.scvi", features = "RGS5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
# FeaturePlot(seu_imm, reduction = "umap.scvi", features = "NOTCH3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "LUM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
# FeaturePlot(seu_imm, reduction = "umap.scvi", features = "DCN", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
# FeaturePlot(seu_imm, reduction = "umap.scvi", features = "PLVAP", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
# FeaturePlot(seu_imm, reduction = "umap.scvi", features = "ACKR1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
# FeaturePlot(seu_imm, reduction = "umap.scvi", features = "ENG", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "PTPRC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "CD68", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "CD3D", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "CD79A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "CD8A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "IL2RA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "KLRD1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "GZMK", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "MS4A1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "TNFRSF17", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "KIT", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "PLD4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "CD163", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm, reduction = "umap.scvi", features = "ICOS", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

seu_imm <- subset(seu_imm, idents = c(36), invert = TRUE)

saveRDS(seu_imm, file = "RDSfiles/seu_imm_DF.RDS")

seu_imm <- JoinLayers(seu_imm)
c18markers <- FindMarkers(seu_imm, ident.1 = 18)

# subset tumor samples, then stratify according to ARID1A expression----
seu_imm_t <- subset(seu_imm, subset = Tumor =="Tumor")
saveRDS(seu_imm_t, file = "RDSfiles/seu_imm_t_DF.RDS")

Idents(seu_imm_t) <- "sample_ID"
seu_imm_t <- subset(seu_imm_t, idents = 15, invert = T) # remove sample_ID 15, which is removed from epithelial analysis (almost no epithelial cells)
Idents(seu_imm_t) <- "sample_ID"
ARID1A_H <- WhichCells(seu_imm_t, ident = c(13,	30,	2,	29,	24,	18,	14,	20,	34,	19,	36,	16,	17,	22))
ARID1A_L <- WhichCells(seu_imm_t, ident = c(38,	33,	3,	8,	27,	28,	12,	5,	40,	39,	7,	32,	26,	10))
seu_imm_t <- SetIdent(seu_imm_t, cells = ARID1A_H, value = "ARID1A_H")
seu_imm_t$ARID1A_level <- Idents(seu_imm_t)
seu_imm_t <- SetIdent(seu_imm_t, cells = ARID1A_L, value = "ARID1A_L")
seu_imm_t$ARID1A_level <- Idents(seu_imm_t)
seu_imm_t$ARID1A_level <- factor(seu_imm_t$ARID1A_level, levels = c("ARID1A_L", "ARID1A_H"))
DimPlot(seu_imm_t, group.by = "ARID1A_level") + NoAxes()

# annodate clusters with Tcells, Bcells and Myeloids----
Idents(seu_imm_t) <- "seurat_clusters"
seu_imm_t <- FindClusters(seu_imm_t, resolution = 0.2, cluster.name = "scvi_clusters", verbose = FALSE)
DimPlot(seu_imm_t, label = TRUE, repel = TRUE) + NoAxes() 
seu_imm_t <- RenameIdents(seu_imm_t, 
                           `0` = "Tcells", `1` = "Tcells", `2` = "Bcells", `3` = "Myeloids", `4` = "Bcells", `5` = "Bcells", 
                           `6` = "Tcells", `7` = "Myeloids", `8` = "Tcells", `9` = "Tcells", `10` = "Tcells", `11` = "Myeloids",`12` = "Tcells")
Idents(seu_imm_t) <- factor(Idents(seu_imm_t), levels(Idents(seu_imm_t))[c(1,2,3)])
DimPlot(seu_imm_t, label = TRUE, repel = TRUE) + NoAxes() 
seu_imm_t$celltype <- Idents(seu_imm_t)

FeaturePlot(seu_imm_t, features = "CD3D", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

saveRDS(seu_imm_t, file = "RDSfiles/seu_imm_t_DF.RDS")

# module scores from Immunity 2023 (not very useful)----
seu_imm_t <- FindClusters(seu_imm_t, resolution = 2, cluster.name = "scvi_clusters", verbose = FALSE)
DimPlot(seu_imm_t, label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('MS4A1', 'CD19', 'VPREB3', 'CD79A', 'BANK1', 'CD79B', 'CD22')),name = "Bcell_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('CLEC9A', 'IDO1', 'CPNE3', 'BATF3')),name = "DC1_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('ID3', 'ENTPD1', 'GZMA', 'CD247', 'CD7', 'HOPX')),name = "IEL_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('ALDOC', 'LINC00299', 'LST1', 'IL4I1', 'AREG')),name = "ILC_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('CD163', 'C1QC', 'C1QA', 'C1QB')),name = "Macrophage_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('TPSAB1', 'CPA3', 'CTSG', 'HDC', 'GATA2', 'VWA5A', 'SLC18A2')),name = "MastCell_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('LAMP3', 'FSCN1', 'CCL19', 'CCL22', 'IDO1', 'CCR7', 'MARCKSL1')),name = "MatureDC_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('IGJ', 'MZB1', 'IGLL5', 'DERL3', 'SSR4', 'TNFRSF17', 'FKBP11', 'SEC11C', 'ANKRD28', 'AL928768.3')),name = "PlasmaCell_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('CTLA4', 'TIGIT', 'TBC1D4', 'BATF', 'TNFRSF4')),name = "Treg_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('CD4', 'FOSB', 'IL7R', 'RORA', 'CD2')),name = "Tcell_CD4_FOSB_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('IL17A', 'IL22', 'CXCR6', 'CCL20')),name = "Tcell_CD4_IL17_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('CD8A', 'CD8B')),name = "Tcell_CD8_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('KLRG1', 'GZMH', 'IFNG', 'CD8B', 'CD8A')),name = "Tcell_CD8_KLRG1_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('CCR7', 'SELL', 'TCF7')),name = "Tcell_Naive_CD4_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('OGT', 'MIAT', 'CELF2', 'RORA', 'ANKRD44', 'ARAP2', 'AKNA', 'CBLB')),name = "Tcell_OGT_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('KLRF1', 'NCAM1', 'KLRD1')),name = "NKcell_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('CHI3L1', 'CYP27A1')),name = "Monocyte_CHI3L1_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('S100A9', 'S100A8', 'FCN1', 'G0S2', 'EREG', 'FPR1')),name = "Monocyte_S100A9_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('CCL3', 'CCL4', 'DAB2', 'A2M')),name = "Macrophage_CCL3_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('CXCL10', 'CXCL9', 'GBP1', 'CXCL11')),name = "Macrophage_CXCL9_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('LYVE1', 'F13A1', 'CCL18')),name = "Macrophage_LYVE1_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('MT1G', 'MT1X', 'MT2A', 'MT1H', 'MT1E', 'MT1F', 'MT1M')),name = "Macrophage_metallothionein_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('PLA2G2D', 'MMP9', 'PTGDS')),name = "Macrophage_PLA2G2D_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('FCER1A', 'CLEC10A')),name = "DC2_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('STMN1', 'HMGB2', 'TCL1A', 'NUSAP1', 'KIAA0101', 'TOP2A', 'TYMS', 'CDK1', 'UBE2C', 'PTTG1')),name = "Cycling_")
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c('CHGA', 'NTS', 'PYY', 'GCG', 'CCK')),name = "Lcell_")

FeaturePlot(seu_imm_t, features = "Bcell_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "DC1_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "IEL_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "ILC_1", cols = c("lightgrey","red","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "Macrophage_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "MastCell_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "MatureDC_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "PlasmaCell_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "Tcell_CD4_FOSB_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "Tcell_CD4_IL17_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "Tcell_CD8_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "Tcell_CD8_KLRG1_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "Tcell_Naive_CD4_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "Tcell_OGT_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "NKcell_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "Monocyte_CHI3L1_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "Monocyte_S100A9_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "Macrophage_CCL3_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "Macrophage_CXCL9_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "Macrophage_LYVE1_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "Macrophage_metallothionein_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "Macrophage_PLA2G2D_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "DC2_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "Cycling_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_imm_t, features = "Lcell_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()


# DA analysis by speckle----
library(speckle)
library(scater)
library(patchwork)

seu_imm_t <- FindClusters(seu_imm_t, resolution = 0.2, cluster.name = "scvi_clusters", verbose = FALSE)
propres <- propeller(clusters = seu_imm_t$seurat_clusters, sample = seu_imm_t$sample_ID, group = seu_imm_t$ARID1A_level)
write.table(propres, file = "results/DA//propeller_imm_t_DF_res0.2.txt", sep ="\t", col.names = T,row.names = F)
plotCellTypeProps(clusters = seu_imm_t$seurat_clusters, sample = seu_imm_t$sample_ID)

# DA analysis by miloR----
library(miloR)
library(SingleCellExperiment)
library(scater)

sce_imm_t <- as.SingleCellExperiment(seu_imm_t)
milo_imm_t <- Milo(sce_imm_t)
milo_imm_t <- buildGraph(milo_imm_t, k = 100, d = 30)
milo_imm_t <- makeNhoods(milo_imm_t, prop = 0.2, k = 100, d=30, refined = TRUE)
plotNhoodSizeHist(milo_imm_t)

milo_imm_t <- countCells(milo_imm_t, meta.data = data.frame(colData(milo_imm_t)), samples="sample_ID")
head(nhoodCounts(milo_imm_t))

design <- data.frame(colData(milo_imm_t))[,c("sample_ID", "ARID1A_level")]
design <- distinct(design)
rownames(design) <- design$sample_ID
design <- design[colnames(nhoodCounts(milo_imm_t)), , drop=FALSE]
design

milo_imm_t <- calcNhoodDistance(milo_imm_t, d=30)
rownames(design) <- design$sample_ID
da_results <- testNhoods(milo_imm_t, design = ~ ARID1A_level, design.df = design)
da_results %>%
  arrange(- SpatialFDR) %>%
  head() 

milo_imm_t <- buildNhoodGraph(milo_imm_t)

plotUMAP(milo_imm_t) + plotNhoodGraphDA(milo_imm_t, da_results, alpha=0.05) +
  plot_layout(guides="collect")

#### manual annotation of clusters----
#### subset T cells, B cells, or Myeloids----
seu_Tcells <- subset(seu_imm_t, subset = celltype == "Tcells")
seu_Bcells <- subset(seu_imm_t, subset = celltype == "Bcells")
seu_Myeloids <- subset(seu_imm_t, subset = celltype == "Myeloids")

#### annotate T cells----
# recluster T cells
seu_Tcells <-DietSeurat(seu_Tcells) 
seu_Tcells <- JoinLayers(seu_Tcells)
seu_Tcells[["RNA"]]$data <- NULL
seu_Tcells[["RNA"]]$scale.data <- NULL
seu_Tcells[["RNA"]] <- split(seu_Tcells[["RNA"]], f = seu_Tcells$sample_ID)
seu_Tcells

seu_Tcells <- NormalizeData(seu_Tcells, verbose = FALSE)
seu_Tcells <- FindVariableFeatures(seu_Tcells, verbose = FALSE)
seu_Tcells <- ScaleData(seu_Tcells, verbose = FALSE)
seu_Tcells <- RunPCA(seu_Tcells, npcs = 30, verbose = FALSE)

seu_Tcells <- IntegrateLayers(
  object = seu_Tcells, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/Users/kinoshitahiroto/miniconda3/envs/scvi-env",
  verbose = TRUE
)

seu_Tcells <- FindNeighbors(seu_Tcells, reduction = "integrated.scvi", dims = 1:30, verbose = FALSE)
seu_Tcells <- FindClusters(seu_Tcells, resolution = 3, cluster.name = "scvi_clusters", verbose = FALSE)
seu_Tcells <- RunUMAP(seu_Tcells, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", verbose = FALSE)

DimPlot(seu_Tcells, label = TRUE, repel = TRUE) + NoAxes()
# DimPlot(seu_Tcells, label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
DimPlot(seu_Tcells, group.by = "sample_ID") + NoAxes()

FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "CD3D", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "CD8A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "CD4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "SELL", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "CD69", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "CXCR3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "IFNG", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "IL4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "IL22", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "RORA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "IL17A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "MS4A1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# NK
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "EOMES", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "PRF1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "NKG7", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "GZMA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# Trrg
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "FOXP3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "CTLA4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "TIGIT", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "IL2RA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# T follicular helper
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "CXCR5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "PDCD1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# CD8 memory
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "FGFBP2", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "S1PR5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "CX3CR1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# MAIT
# FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "TRAV1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "TRAV2", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "SLC4A10", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# gdT
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "TRGV2", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "TRGV4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "TRGV5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "TRGV7", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# ILC3 etc.
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "RORC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "IL1R1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "IL23R", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "KIT", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "TNFSF4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "PCDH9", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "TNFRSF11A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "TNFSF11", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "NCR1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "NCR2", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "ITGAE", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "CCR9", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# ILC2
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "PTGDR2", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "HPGDS", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "IL1RL1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Tcells, reduction = "umap.scvi", features = "KRT1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# annotate Tcell clusters
Idents(seu_Tcells) <- "seurat_clusters"
seu_Tcells <- subset(seu_Tcells, idents = 44, invert = T) # cluster#44 looks like contamination of B cells
seu_Tcells <- RenameIdents(seu_Tcells, 
                          `0` = "CD8", `1` = "Treg", `2` = "CD8", `3` = "CD8", `4` = "CD8", `5` = "CD8", `6` = "Tmem", `7` = "CD4", `8` = "CD4", `9` = "CD4", 
                          `10` = "CD8", `11` = "CD8",`12` = "NK", `13` = "CD8", `14` = "CD4",`15` = "CD8", `16` = "CD4", `17` = "MAIT",`18` = "CD4", `19` = "CD8",
                          `20` = "CD8", `21` = "CD8",`22` = "CD4", `23` = "CD8", `24` = "NK",`25` = "Treg", `26` = "TfH", `27` = "CD8",`28` = "CD8", `29` = "Tmem",
                          `30` = "Treg", `31` = "NKT",`32` = "CD4+CD8+", `33` = "NK", `34` = "Treg",`35` = "Treg", `36` = "CD8", `37` = "CD4-CD8-",`38` = "NK", `39` = "CD8",
                          `40` = "CD8", `41` = "CD8",`42` = "CD8", `43` = "NKT", `45` = "NK", `46` = "CD4", `47` = "Tmem",`48` = "CD4-CD8-", `49` = "CD4")
levels(Idents(seu_Tcells))
Idents(seu_Tcells) <- factor(Idents(seu_Tcells), levels(Idents(seu_Tcells))[c(10,4,1,9,2,5,8,3,6,7)])
DimPlot(seu_Tcells, label = TRUE, repel = TRUE) + NoAxes() 
seu_Tcells$cellSubtype <- Idents(seu_Tcells)
saveRDS(seu_Tcells, file = "RDSfiles/seu_Tcells_annotated_DF.RDS")

#### annotate B cells----
# recluster B cells
seu_Bcells <-DietSeurat(seu_Bcells) 
seu_Bcells <- JoinLayers(seu_Bcells)
seu_Bcells[["RNA"]]$data <- NULL
seu_Bcells[["RNA"]]$scale.data <- NULL
seu_Bcells[["RNA"]] <- split(seu_Bcells[["RNA"]], f = seu_Bcells$sample_ID)
seu_Bcells

seu_Bcells <- NormalizeData(seu_Bcells, verbose = FALSE)
seu_Bcells <- FindVariableFeatures(seu_Bcells, verbose = FALSE)
seu_Bcells <- ScaleData(seu_Bcells, verbose = FALSE)
seu_Bcells <- RunPCA(seu_Bcells, npcs = 30, verbose = FALSE)

seu_Bcells <- IntegrateLayers(
  object = seu_Bcells, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/Users/kinoshitahiroto/miniconda3/envs/scvi-env",
  verbose = TRUE
)

seu_Bcells <- FindNeighbors(seu_Bcells, reduction = "integrated.scvi", dims = 1:30, verbose = FALSE)
seu_Bcells <- FindClusters(seu_Bcells, resolution = 3, cluster.name = "scvi_clusters", verbose = FALSE)
seu_Bcells <- RunUMAP(seu_Bcells, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", verbose = FALSE)

DimPlot(seu_Bcells, label = TRUE, repel = TRUE) + NoAxes()
# DimPlot(seu_Bcells, label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
DimPlot(seu_Bcells, group.by = "sample_ID") + NoAxes()

FeaturePlot(seu_Bcells, reduction = "umap.scvi", features = "MS4A1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Bcells, reduction = "umap.scvi", features = "CD19", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Bcells, reduction = "umap.scvi", features = "CD38", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Bcells, reduction = "umap.scvi", features = "IGHA1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Bcells, reduction = "umap.scvi", features = "IGHG1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Bcells, reduction = "umap.scvi", features = "IGHM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Bcells, reduction = "umap.scvi", features = "FCMR", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Bcells, reduction = "umap.scvi", features = "SELL", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Bcells, reduction = "umap.scvi", features = "FCRL4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Bcells, reduction = "umap.scvi", features = "SDC1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Bcells, reduction = "umap.scvi", features = "MKI67", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Bcells, reduction = "umap.scvi", features = "HMGB2", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Bcells, reduction = "umap.scvi", features = "TUBA1B", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Bcells, reduction = "umap.scvi", features = "TNFRSF17", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# annotate Bcell clusters
Idents(seu_Bcells) <- "seurat_clusters"
seu_Bcells <- RenameIdents(seu_Bcells, 
                           `0` = "PlasmaIgA", `1` = "B1", `2` = "PlasmaIgA", `3` = "PlasmaIgA", `4` = "B1", `5` = "PlasmaIgA", `6` = "PlasmaIgA", `7` = "PlasmaIgA", `8` = "PlasmaIgA", `9` = "PlasmaIgA", 
                           `10` = "B1", `11` = "PlasmaIgA",`12` = "NaiveB", `13` = "PlasmaIgA", `14` = "PlasmaIgG",`15` = "PlasmaIgA", `16` = "PlasmaIgA", `17` = "PlasmaIgA",`18` = "PlasmaIgA", `19` = "PlasmaIgA",
                           `20` = "PlasmaIgG", `21` = "PlasmaIgA",`22` = "PlasmaIgA", `23` = "PlasmaIgA", `24` = "PlasmaIgG",`25` = "PlasmaIgA", `26` = "PlasmaIgA", `27` = "PlasmaIgG",`28` = "B2", `29` = "PlasmaIgA",
                           `30` = "PlasmaIgA", `31` = "B1",`32` = "PlasmaIgA", `33` = "CyclingB", `34` = "Bmem",`35` = "B1", `36` = "PlasmaIgA")
levels(Idents(seu_Bcells))
Idents(seu_Bcells) <- factor(Idents(seu_Bcells), levels(Idents(seu_Bcells))[c(3,2,5,6,7,1,4)])
DimPlot(seu_Bcells, label = TRUE, repel = TRUE) + NoAxes() 
seu_Bcells$cellSubtype <- Idents(seu_Bcells)
saveRDS(seu_Bcells, file = "RDSfiles/seu_Bcells_annotated_DF.RDS")


#### annotate Myeloids----
# recluster Myeloids
seu_Myeloids <-DietSeurat(seu_Myeloids) 
seu_Myeloids <- JoinLayers(seu_Myeloids)
seu_Myeloids[["RNA"]]$data <- NULL
seu_Myeloids[["RNA"]]$scale.data <- NULL
seu_Myeloids[["RNA"]] <- split(seu_Myeloids[["RNA"]], f = seu_Myeloids$sample_ID)
seu_Myeloids

seu_Myeloids <- NormalizeData(seu_Myeloids, verbose = FALSE)
seu_Myeloids <- FindVariableFeatures(seu_Myeloids, verbose = FALSE)
seu_Myeloids <- ScaleData(seu_Myeloids, verbose = FALSE)
seu_Myeloids <- RunPCA(seu_Myeloids, npcs = 30, verbose = FALSE)

seu_Myeloids <- IntegrateLayers(
  object = seu_Myeloids, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/Users/kinoshitahiroto/miniconda3/envs/scvi-env",
  verbose = TRUE
)

seu_Myeloids <- FindNeighbors(seu_Myeloids, reduction = "integrated.scvi", dims = 1:30, verbose = FALSE)
seu_Myeloids <- FindClusters(seu_Myeloids, resolution = 3, cluster.name = "scvi_clusters", verbose = FALSE)
seu_Myeloids <- RunUMAP(seu_Myeloids, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", verbose = FALSE)

DimPlot(seu_Myeloids, label = TRUE, repel = TRUE) + NoAxes()
# DimPlot(seu_Myeloids, label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
DimPlot(seu_Myeloids, group.by = "sample_ID") + NoAxes()

FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "CD68", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "FCN1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "S100A4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "S100A6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "GATA2", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "CPA3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "HPGDS", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "CLEC9A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "CLEC10A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "LAMP3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "PLD4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "CLEC4C", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "JCHAIN", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "CD163", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "C1QB", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "C1QC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "LYVE1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "RNASE1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "SPP1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "MMP9", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_Myeloids, reduction = "umap.scvi", features = "GPNMB", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# annotate Myeloid clusters
Idents(seu_Myeloids) <- "seurat_clusters"
seu_Myeloids <- subset(seu_Myeloids, idents = 27, invert = T) # clusters #27 looks like contamination of plasma cells 
seu_Myeloids <- RenameIdents(seu_Myeloids, 
                           `0` = "Mast", `1` = "cDC2", `2` = "Monocyte", `3` = "infl.Macrophage", `4` = "Mast", `5` = "Monocyte", `6` = "Mast", `7` = "infl.Macrophage", `8` = "Macrophage", `9` = "Macrophage", 
                           `10` = "infl.Macrophage", `11` = "infl.Macrophage",`12` = "Monocyte", `13` = "Mast", `14` = "Macrophage",`15` = "lymphoidDC", `16` = "Macrophage", `17` = "Monocyte",`18` = "Macrophage", `19` = "Macrophage",
                           `20` = "Monocyte", `21` = "Mast",`22` = "Macrophage", `23` = "cDC1", `24` = "Macrophage",`25` = "LYVE1Macrophage", `26` = "Mast",`28` = "pDC", `29` = "Macrophage",
                           `30` = "Monocyte", `31` = "DC?",`32` = "infl.Macrophage", `33` = "Monocyte")
levels(Idents(seu_Myeloids))
Idents(seu_Myeloids) <- factor(Idents(seu_Myeloids), levels(Idents(seu_Myeloids))[c(1,3,5,4,8,7,2,6,9,10)])
DimPlot(seu_Myeloids, label = TRUE, repel = TRUE) + NoAxes() 
seu_Myeloids$cellSubtype <- Idents(seu_Myeloids)
saveRDS(seu_Myeloids, file = "RDSfiles/seu_Myeloids_annotated_DF.RDS")

#### combine annotated seurat objects and do DA analysis----
seu_Tcells <- JoinLayers(seu_Tcells)
seu_Bcells <- JoinLayers(seu_Bcells)
seu_Myeloids <- JoinLayers(seu_Myeloids)
seu_imm_t <- merge(x = seu_Tcells, y = list(seu_Bcells, seu_Myeloids))

seu_imm_t <-DietSeurat(seu_imm_t) 
seu_imm_t <- JoinLayers(seu_imm_t)
seu_imm_t[["RNA"]]$data <- NULL
seu_imm_t[["RNA"]]$scale.data <- NULL
seu_imm_t[["RNA"]] <- split(seu_imm_t[["RNA"]], f = seu_imm_t$sample_ID)
seu_imm_t

seu_imm_t <- NormalizeData(seu_imm_t, verbose = FALSE)
seu_imm_t <- FindVariableFeatures(seu_imm_t, verbose = FALSE)
seu_imm_t <- ScaleData(seu_imm_t, verbose = FALSE)
seu_imm_t <- RunPCA(seu_imm_t, npcs = 30, verbose = FALSE)

seu_imm_t <- IntegrateLayers(
  object = seu_imm_t, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/Users/kinoshitahiroto/miniconda3/envs/scvi-env",
  verbose = TRUE
)

seu_imm_t <- FindNeighbors(seu_imm_t, reduction = "integrated.scvi", dims = 1:30, verbose = FALSE)
seu_imm_t <- FindClusters(seu_imm_t, resolution = 2, cluster.name = "scvi_clusters", verbose = FALSE)
seu_imm_t <- RunUMAP(seu_imm_t, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", verbose = FALSE)

Idents(seu_imm_t) <- "cellSubtype"
DimPlot(seu_imm_t, label = TRUE, repel = TRUE) + NoAxes()
saveRDS(seu_imm_t, file = "RDSfiles/seu_imm_t_DF_annotated.RDS")

#### additional analysis----
# ILC2 signature from nat commun 2023 by OKeefe et al.----
seu_imm_t <- JoinLayers(seu_imm_t)
seu_imm_t <- AddModuleScore(seu_imm_t, features = list(c("GATA3", "IL13", "ICOS", "KLRG1", "CRTH2", "IL5", "IL4")),name = "ILC2_")
VlnPlot(seu_imm_t, features = "ILC2_1", group.by = "ARID1A_level", pt.size = 0.0)+ 
  geom_signif(comparisons = list(c("ARID1A_H","ARID1A_L")),map_signif_level = function(x) paste("p =", scales::pvalue(x)), textsize =4) + 
  ylim(NA,0.1) +
  stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(seu_imm_t, features = "ILC2_1", group.by = "celltype", split.by = "ARID1A_level", pt.size = 0.0)
Idents(seu_imm_t) <- "cellSubtype"
FeaturePlot(seu_imm_t, features = "ILC2_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# pseudobulk like analysis manually----
Idents(seu_imm_t) <- "sample_ID"

# make data.frame of sample_ID, ARID1A_level, mean of ILC2 signature
## vector of mean of ILC2 signature with sample_ID
mean_sig <-c()
for (id in  as.integer(levels(factor(seu_imm_t$sample_ID)))) {
  mean_sig <- c(mean_sig, mean(seu_imm_t$ILC2_1[which(colnames(seu_imm_t) %in%  WhichCells(seu_imm_t, idents = id))]))
}
sample_ID <- as.integer(levels(factor(seu_imm_t$sample_ID)))
df <- cbind(sample_ID, mean_sig) %>% as.data.frame()

## ARID1A_level H or L by sample_ID (from average expression in epithelial cells, refer to "Kumar_scvi_epi.R")
ARID1A_H <- c(13,	30,	2,	29,	24,	18,	14,	20,	34,	19,	36,	16,	17,	22)
ARID1A_L <- c(38,	33,	3,	8,	27,	28,	12,	5,	40,	39,	7,	32,	26,	10)
sample_ID <- c(ARID1A_H, ARID1A_L)
ARID1A_level <- c(rep("H", 14), rep("L", 14))
ARID1A_level <- cbind(sample_ID, ARID1A_level) %>% as.data.frame()
ARID1A_level$sample_ID <- ARID1A_level$sample_ID %>% as.double()

## join two data.frame
df <- left_join(df, ARID1A_level, by = "sample_ID")

## violin plot
ggplot(df, aes(x = ARID1A_level, y = mean_sig, fill = ARID1A_level)) + 
  geom_violin() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  geom_signif(comparisons = list(c("H","L")),map_signif_level = function(x) paste("p =", scales::pvalue(x)), textsize =4) + 
  ylim(NA,0.06) + 
  stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95) +
  theme_classic()

# the same analysis as above in Tcells only----
seu_Tcells <- JoinLayers(seu_Tcells)
seu_Tcells <- AddModuleScore(seu_Tcells, features = list(c("GATA3", "IL13", "ICOS", "KLRG1", "CRTH2", "IL5", "IL4")),name = "ILC2_")
VlnPlot(seu_Tcells, features = "ILC2_1", group.by = "ARID1A_level", pt.size = 0.0) + 
  geom_signif(comparisons = list(c("ARID1A_H","ARID1A_L")),map_signif_level = function(x) paste("p =", scales::pvalue(x)), textsize =4) + 
  ylim(NA,0.1) +
  stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95)
VlnPlot(seu_Tcells, features = "ILC2_1", group.by = "cellSubtype", split.by = "ARID1A_level", pt.size = 0.0)
Idents(seu_Tcells) <- "cellSubtype"
FeaturePlot(seu_Tcells, features = "ILC2_1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
Idents(seu_Tcells) <- "sample_ID"

# make data.frame of sample_ID, ARID1A_level, mean of ILC2 signature
## vector of mean of ILC2 signature with sample_ID
Idents(seu_Tcells) <- "sample_ID"
sample_ID <- as.integer(levels(factor(seu_Tcells$sample_ID)))
mean_sig <-c()
for (id in sample_ID) {
  mean_sig <- c(mean_sig, mean(seu_Tcells$ILC2_1[which(colnames(seu_Tcells) %in%  WhichCells(seu_Tcells, idents = id))]))
}

df <- cbind(sample_ID, mean_sig) %>% as.data.frame()

## ARID1A_level H or L by sample_ID (from average expression in epithelial cells, refer to "Kumar_scvi_epi.R")
ARID1A_H <- c(13,	30,	2,	29,	24,	18,	14,	20,	34,	19,	36,	16,	17,	22)
ARID1A_L <- c(38,	33,	3,	8,	27,	28,	12,	5,	40,	39,	7,	32,	26,	10)
sample_ID <- c(ARID1A_H, ARID1A_L)
ARID1A_level <- c(rep("H", 14), rep("L", 14))
ARID1A_level <- cbind(sample_ID, ARID1A_level) %>% as.data.frame()
ARID1A_level$sample_ID <- ARID1A_level$sample_ID %>% as.double()

## join two data.frame
df <- left_join(df, ARID1A_level, by = "sample_ID")

## violin plot
ggplot(df, aes(x = ARID1A_level, y = mean_sig, fill = ARID1A_level)) + 
  geom_violin() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  geom_signif(comparisons = list(c("H","L")),map_signif_level = function(x) paste("p =", scales::pvalue(x)), textsize =4) + 
  ylim(NA,0.03) + 
  stat_summary(fun.y=median, geom="point", size=5, colour="black", shape=95) +
  theme_classic()
