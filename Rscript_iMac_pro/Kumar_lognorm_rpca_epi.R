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

seu_epi <- readRDS("RDSfiles/seu_lognorm_rpca_epi.RDS")


seu_epi <- DietSeurat(seu_epi) 
seu_epi <- JoinLayers(seu_epi)
seu_epi[["RNA"]]$data <- NULL
seu_epi[["RNA"]]$scale.data <- NULL
seu_epi[["RNA"]] <- split(seu_epi[["RNA"]], f = seu_epi$sample_ID)
seu_epi #sample 31 is missing (probably there is no epithelial cells in this samples)

seu_epi <- NormalizeData(seu_epi, verbose = FALSE)
seu_epi <- FindVariableFeatures(seu_epi, verbose = FALSE)
seu_epi <- ScaleData(seu_epi, verbose = FALSE)
seu_epi <- RunPCA(seu_epi, npcs = 30, verbose = FALSE)

table(seu_epi$sample_ID) # check cell numbers in each sample
# turned out sample31 was not included (no epithelial cells in sample31)
# choose samples with >100 epithelial cells
ncell_df <- as.data.frame(table(seu_epi$sample_ID))
ncell100_df <- filter(ncell_df, Freq > 100)
ncell100_ID <- ncell100_df[,1] %>% as.character %>% as.double
seu_epi <- subset(seu_epi, subset = sample_ID %in% ncell100_ID)


seu_epi <- IntegrateLayers(
  object = seu_epi, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.rpca",
  reference = c(4, 5, 7, 23),
  k.weight = 90,
  verbose = FALSE
)

seu_epi <- FindNeighbors(seu_epi, reduction = "integrated.rpca", dims = 1:30, verbose = FALSE)
seu_epi <- FindClusters(seu_epi, resolution = 1, cluster.name = "spca_clusters", verbose = FALSE)
seu_epi <- RunUMAP(seu_epi, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca", verbose = FALSE)

DimPlot(seu_epi, label = TRUE, repel = TRUE) + NoLegend() + NoAxes()
DimPlot(seu_epi, group.by = "sample_ID") + NoAxes()
DimPlot(seu_epi, group.by = "Tumor") + NoAxes()

# saveRDS(seu_epi, file = "RDSfiles/seu_lognorm_rpca_epi.RDS")

FeaturePlot(seu_epi,features = "EPCAM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "PTPRC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "FN1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "MUC5AC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "TFF1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "LIPF", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "PGA3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "REG4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "TFF3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "CHGA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "RGS5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "NOTCH3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "LUM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "DCN", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "PLVAP", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "ACKR1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "CD79A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "CD8A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "IL2RA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "KLRD1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "MS4A1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "TNFRSF17", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "KIT", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu_epi,features = "LRG1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "CD38", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "ENG", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "MUC6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "ATP4B", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

seu_epi <- subset(seu_epi, idents = c(0:11,13,15,16,18,20,21))

# repeat clustering and check feature plots until 

saveRDS(seu_epi, file = "RDSfiles/seu_lognorm_rpca_epi_2.RDS")

VlnPlot(seu_epi, features = "LRG1", group.by = "seurat_clusters", split.by = "Tumor", pt.size = 0)

