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
  conda_env = "/usr/local/Caskroom/miniforge/base/envs/scvi-env2",
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

Idents(seu_imm_t) <- "sample_ID"
seu_imm_t <- subset(seu_imm_t, idents = 15, invert = T) # remove sample_ID 15, which is removed from epithelial analysis (almost no epithelial cells)
ARID1A_H <- WhichCells(seu_imm_t, ident = c(13,	30,	2,	29,	24,	18,	14,	20,	34,	19,	36,	16,	17,	22))
ARID1A_L <- WhichCells(seu_imm_t, ident = c(38,	33,	3,	8,	27,	28,	12,	5,	40,	39,	7,	32,	26,	10))
seu_imm_t <- SetIdent(seu_imm_t, cells = ARID1A_H, value = "ARID1A_H")
seu_imm_t$ARID1A_level <- Idents(seu_imm_t)
seu_imm_t <- SetIdent(seu_imm_t, cells = ARID1A_L, value = "ARID1A_L")
seu_imm_t$ARID1A_level <- Idents(seu_imm_t)
seu_imm_t$ARID1A_level <- factor(seu_imm_t$ARID1A_level, levels = c("ARID1A_L", "ARID1A_H"))
DimPlot(seu_imm_t, group.by = "ARID1A_level") + NoAxes()

saveRDS(seu_imm_t, file = "RDSfiles/seu_imm_t_DF.RDS")

