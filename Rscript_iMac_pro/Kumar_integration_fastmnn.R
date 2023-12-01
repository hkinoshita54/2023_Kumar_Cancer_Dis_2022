library(tidyverse)
library(cowplot)
library(ggplot2)
library(patchwork)
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratWrappers)
library(reticulate)
options(future.globals.maxSize = 1e10)
library(harmony)


# the standard (lognormalization) pipeline from the vignette
seu <- readRDS(file = "RDSfiles/seu_DF.RDS")
seu <- JoinLayers(seu)
# seu <- DietSeurat(seu)
seu[["RNA"]]$data <- NULL
seu[["RNA"]]$scale.data <- NULL
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$sample_ID)
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)

# Dimplot with integration
seu <- IntegrateLayers(
  object = seu, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)
seu <- FindNeighbors(seu, reduction = "integrated.mnn", dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 1, cluster.name = "mnn_clusters", verbose = FALSE)
seu <- RunUMAP(seu, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn", verbose = FALSE)
DimPlot(seu, reduction = "umap.mnn", group.by = c("mnn_clusters"), label = TRUE, repel = TRUE) + NoLegend() + NoAxes()
# DimPlot(seu, reduction = "umap.scvi", group.by = c("sample_ID", "scvi_clusters"), combine = FALSE)
Idents(seu)
DimPlot(seu, reduction = "umap.mnn", group.by = "sample_ID") + NoAxes()

saveRDS(seu, file = "RDSfiles/seu_DF.RDS")

# Check which cluster is epithelial, immune or stroma ----

FeaturePlot(seu, reduction = "umap.mnn", features = "EPCAM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "PTPRC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "FN1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu, reduction = "umap.mnn", features = "MUC5AC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "TFF1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "MUC6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "LIPF", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "PGA3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "REG4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "TFF3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "CHGA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu, reduction = "umap.mnn", features = "RGS5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "NOTCH3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "LUM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "DCN", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "PLVAP", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "ACKR1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu, reduction = "umap.mnn", features = "CD68", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "CD3D", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "CD79A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu, reduction = "umap.mnn", features = "IL2RA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "KLRD1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "MS4A1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "TNFRSF17", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "KIT", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu, reduction = "umap.mnn", features = "LRG1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "CD38", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu, reduction = "umap.mnn", features = "ENG", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu, reduction = "umap.mnn", features = "TRPM5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# subset epithelial, stromal and immune clusters
seu_epi <- subset(seu, idents = c(7,13,17,20,23,28)) # c18 and c37 positive fro both immune and epithelial markers, mostly from sample 29 > omit
# seu_str <- subset(seu, idents = c(6,20,23,26,33,36,41,48,49))
# seu_imm <- subset(seu, idents = c(0:5,7:9,12:14,17,21,22,25,27,30,35,39,40,42:44,46,47,51)) # c45 and c52 negative for all the markers > omit

# saveRDS(seu_epi, file = "RDSfiles/seu_epi_DF.RDS")
# saveRDS(seu_str, file = "RDSfiles/seu_str_DF.RDS")
# saveRDS(seu_imm, file = "RDSfiles/seu_imm_DF.RDS")