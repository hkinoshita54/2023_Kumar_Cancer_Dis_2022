library(tidyverse)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratWrappers)
library(reticulate)
options(future.globals.maxSize = 1e10)
library(msigdbr)
library(fgsea)
library(DESeq2)

# re-integrate epithelial subset and clean it to remove immune or stroma contaminatino----
# latest analysis if from Kumar_DF, in which mt<12 cells were selected
seu_epi <- readRDS("RDSfiles/seu_epi_scvi.RDS")

seu_epi <-DietSeurat(seu_epi) 
seu_epi <- JoinLayers(seu_epi)
seu_epi[["RNA"]]$data <- NULL
seu_epi[["RNA"]]$scale.data <- NULL
seu_epi[["RNA"]] <- split(seu_epi[["RNA"]], f = seu_epi$sample_ID)
seu_epi #sample 11 and 31 are missing (probably there are no epithelial cells in those samples)

seu_epi <- NormalizeData(seu_epi, verbose = FALSE)
seu_epi <- FindVariableFeatures(seu_epi, verbose = FALSE)
seu_epi <- ScaleData(seu_epi, verbose = FALSE)
seu_epi <- RunPCA(seu_epi, npcs = 30, verbose = FALSE)

seu_epi <- IntegrateLayers(
  object = seu_epi, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)

seu_epi <- FindNeighbors(seu_epi, reduction = "integrated.mnn", dims = 1:30, verbose = FALSE)
seu_epi <- FindClusters(seu_epi, resolution = 1.2, cluster.name = "mnn_clusters", verbose = FALSE)
seu_epi <- RunUMAP(seu_epi, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn", verbose = FALSE)

DimPlot(seu_epi, label = T, repel = TRUE) + NoAxes()
DimPlot(seu_epi, group.by = "sample_ID") + NoAxes()
DimPlot(seu_epi, group.by = "Tumor") + NoAxes()

FeaturePlot(seu_epi,features = "EPCAM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "PTPRC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "VIM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu_epi,features = "MUC5AC", cols = c("lightgrey","darkred"), label = T, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "MUC6", cols = c("lightgrey","darkred"), label = T, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "LIPF", cols = c("lightgrey","darkred"), label = T, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "ATP4B", cols = c("lightgrey","darkred"), label = T, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "TFF3", cols = c("lightgrey","darkred"), label = T, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "CHGA", cols = c("lightgrey","darkred"), label = T, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "RGS5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "DCN", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "PLVAP", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "CD79A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "CD3D", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "CD68", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "TNFRSF17", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "GHRL", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# remove clusters contaminated with stroma or immune
seu_epi <- subset(seu_epi, idents = c(12,14), invert = TRUE)
# repeat clustering and check feature plots until it is clean

saveRDS(seu_epi, file = "RDSfiles/seu_epi_DF.RDS")


