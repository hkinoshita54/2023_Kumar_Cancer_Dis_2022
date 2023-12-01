library(tidyverse)
library(data.table)
library(Matrix)
library(cowplot)
library(ggplot2)
library(patchwork)
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratWrappers)

patients_data <- read.delim("meta_data/patients_data.txt", header = TRUE)
seu_epi <- readRDS(file = "RDSfiles/seu_epi_v5.RDS")
seu_epi <- DietSeurat(seu_epi)

# the standard (lognormalization) pipeline from the vignette
seu_epi <- NormalizeData(seu_epi, verbose = FALSE)
seu_epi <- FindVariableFeatures(seu_epi, verbose = FALSE)
seu_epi <- ScaleData(seu_epi, verbose = FALSE)
seu_epi <- RunPCA(seu_epi, npcs = 30, verbose = FALSE)

# Dimplot without integration
seu_epi <- FindNeighbors(seu_epi, reduction = "pca", dims = 1:30, verbose = FALSE)
seu_epi <- FindClusters(seu_epi, resolution = 0.8, cluster.name = "unintegrated_clusters", verbose = FALSE)
seu_epi <- RunUMAP(seu_epi, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated", verbose = FALSE)
DimPlot(seu_epi, reduction = "umap.unintegrated", group.by = c("sample_ID", "unintegrated_clusters"), combine = FALSE)

# Dimplot with integration (RPCA with reference)
options(future.globals.maxSize = 8000 * 1024^2) # to avoid error related to "future memory"

table(seu_epi$sample_ID) # check cell numbers in each sample
# turned out sample31 was not included (no epithelial cells in sample31)
# choose samples with >30 epithelial cells
ncell_df <- as.data.frame(table(seu_epi$sample_ID))
over30_df <- filter(ncell_df, Freq > 30)
over30_ID <- over30_df[,1] %>% as.character %>% as.double
seu_epi2 <- subset(seu_epi, subset = sample_ID %in% over30_ID)

# choose reference for integration (intestinal tumor, diffuse tumor, intestinal normal, diffuse(mixed) normal)
ncell_vec <- ncell_df[,2]
ncell_vec_full <- c(ncell_vec[1:30], 0, ncell_vec[31:39])
patients_data$epi_cells <- ncell_vec_full
# sample 1, 6, 28, 29

seu_epi2 <- IntegrateLayers(
  object = seu_epi2, method = RPCAIntegration, 
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  # normalization.method = "SCT",
  reference = c(1, 4, 24, 25),
  k.weight = 50,
  verbose = FALSE
)
seu_epi2 <- FindNeighbors(seu_epi2, reduction = "integrated.rpca", dims = 1:30, verbose = FALSE)
seu_epi2 <- FindClusters(seu_epi2, resolution = 1, cluster.name = "rpca_clusters", verbose = FALSE)
seu_epi2 <- RunUMAP(seu_epi2, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca", verbose = FALSE)
DimPlot(seu_epi2, reduction = "umap.rpca",
        group.by = c("sample_ID", 
                     # "unintegrated_clusters", 
                     "rpca_clusters"),
        combine = FALSE
)
# DimPlot(seu_epi, label = TRUE, repel = TRUE) + NoLegend() + NoAxes()
# DimPlot(seu_epi, group.by = "sample_ID") + NoAxes()

saveRDS(seu_epi2, file = "RDSfiles/seu_epi_v5_int.RDS")

FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "EPCAM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "MUC5AC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "TFF1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "LIPF", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "PGA3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "REG4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "TFF3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "CHGA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "TRPM5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "FN1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "RGS5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "NOTCH3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "LUM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "DCN", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "PLVAP", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "ACKR1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "PTPRC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "CD68", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "CD3D", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "CD79A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "CD8A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "IL2RA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "KLRD1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "GZMK", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "MS4A1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "TNFRSF17", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "KIT", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "PLD4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "CD163", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "LRG1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "CD38", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "ENG", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "MUC6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "CD44", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "CLU", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi2, reduction = "umap.rpca", features = "S100A4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()