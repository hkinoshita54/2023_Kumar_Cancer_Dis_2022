library(Seurat)
library(ggplot2)

seu_epi <- readRDS(file = "RDSfiles/seu_epi.RDS")
seu_epi_list <- SplitObject(seu_epi, split.by = "sample_ID")
seu_epi_list <- lapply(X = seu_epi_list, FUN = function(x) {
  x <- DietSeurat(x)
})



seu_epi_2 <- readRDS(file = "RDSfiles/seu_epi_2.RDS")
seu_epi_2_list <- SplitObject(seu_epi, split.by = "sample_ID")
seu_epi_2_list <- lapply(X = seu_epi_2_list, FUN = function(x) {
  x <- DietSeurat(x)
})

# rm(seu_epi, seu_epi_2)

seu_epi_12 <- 

  
  seu_epi_12 <- lapply(X = seu_epi_12, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = seu_epi_12)
seu_epi_12 <- lapply(X = seu_epi_12, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = seu_epi_12, reference = c(1, 3, 4, 12), 
                                  reduction = "rpca", dims = 1:50)
seu_epi_12 <- IntegrateData(anchorset = anchors, dims = 1:50)
seu_epi_12 <- ScaleData(seu_epi_12, verbose = FALSE)
seu_epi_12 <- RunPCA(seu_epi_12, verbose = FALSE)

saveRDS(seu_epi_12, file = "./seu_epi_12.rds")

seu_epi_12 <- RunUMAP(seu_epi_12, dims = 1:30)
seu_epi_12 <- FindNeighbors(seu_epi_12, reduction = "pca", dims = 1:30)
seu_epi_12 <- FindClusters(seu_epi_12, resolution = 0.5)
DimPlot(seu_epi_12, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend() + NoAxes()

