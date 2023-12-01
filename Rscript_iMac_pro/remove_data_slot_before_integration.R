seu_copy <- seu

# test if it is necessary to do Normalize data again after removing data slot
# go directly to FindCariableFeatures after spliitting
seu <- subset(seu_copy, subset = sample_ID %in% c(1:10))
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$sample_ID)
seu <- FindVariableFeatures(seu, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- IntegrateLayers(
  object = seu, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/usr/local/Caskroom/miniforge/base/envs/scvi-env2",
  verbose = TRUE
)
seu <- FindNeighbors(seu, reduction = "integrated.scvi", dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 1, cluster.name = "scvi_clusters", verbose = FALSE)
seu <- RunUMAP(seu, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE) + NoLegend() + NoAxes()

# remove data slot first and then start from Normalize data
seu <- subset(seu_copy, subset = sample_ID %in% c(1:10))
seu[["RNA"]]$data <- NULL
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$sample_ID)
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
seu <- IntegrateLayers(
  object = seu, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/usr/local/Caskroom/miniforge/base/envs/scvi-env2",
  verbose = TRUE
)
seu <- FindNeighbors(seu, reduction = "integrated.scvi", dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 1, cluster.name = "scvi_clusters", verbose = FALSE)
seu <- RunUMAP(seu, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE) + NoLegend() + NoAxes()
# There was a difference from previous dimplot. Probably the second way is safer (remove data slot and do NormalizeData before integration)



seu <- IntegrateLayers(
  object = seu, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)
seu <- FindNeighbors(seu, reduction = "integrated.mnn", dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 1, cluster.name = "mnn_clusters", verbose = FALSE)
seu <- RunUMAP(seu, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn", verbose = FALSE)
DimPlot(seu, reduction = "umap.mnn", group.by = c("mnn_clusters"), label = TRUE, repel = TRUE) + NoLegend() + NoAxes()

seu <- seu_copy
rm(seu_copy)
