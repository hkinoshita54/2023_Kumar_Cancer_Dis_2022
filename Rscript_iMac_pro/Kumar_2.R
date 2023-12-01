library(tidyverse)
library(data.table)
library(cowplot)
library(Seurat)
library(sctransform)
library(ggplot2)
library(Matrix)
library(glmGamPoi)

# make a mata_data table for each saample----
## patients_data.txt was made with excel from suppl. table1 and series_matrix
patients_data <- read.delim("meta_data/patients_data.txt", header = TRUE)
files_list <- list.files(path = "raw_data/GSE183904_RAW_2", pattern = "*.csv.gz", full.names = T)
mtx_list <- list()
for (i in seq_along(files_list)) {
  mtx <- fread(files_list[[i]], header = TRUE)
  mtx <- data.frame(mtx, row.names = 1)
  mtx <- Matrix(as.matrix(mtx), sparse = TRUE)
  mtx_list[[i]] <- mtx
}

meta_list <- list()
for (i in seq_along(mtx_list)) {
  meta <- data.frame(cell_ID = colnames(mtx_list[[i]]), sample_ID = 20 + i)
  meta <- left_join(meta, patients_data)
  meta <- column_to_rownames(meta, var = "cell_ID")
  meta_list[[i]] <- meta
}

# Create Seurat object----
seu_list <- list()
for (i in seq_along(mtx_list)) {
  seu <- CreateSeuratObject(mtx_list[[i]], min.cells = 3, meta.data = meta_list[[i]])
  seu_list[[i]] <- seu
}

# QC for percent.mt and nFeature_RNA----
seu_merge <- merge(seu_list[[1]], seu_list[2:20], add.cell.ids = c(1:20))
seu_merge[["percent.mt"]] <- PercentageFeatureSet(seu_merge, pattern = "^MT-")
VlnPlot(seu_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
## assume cell are already filtered with nFeature_RNA<6000 and percent.mt<20

rm(seu_merge)
gc()

# modified "integrating large datasets" and "introduction to scRNA-seq integration" from Seurat vignettes----

seu_list <- lapply(X = seu_list, FUN = function(x){
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = c("percent.mt"), vst.flavor = "v2", verbose = TRUE)
})

features <- SelectIntegrationFeatures(object.list = seu_list, nfeatures = 3000)
seu_list <- PrepSCTIntegration(object.list = seu_list, anchor.features = features)

rm(meta, meta_list, mtx, mtx_list, seu)
gc()

anchors <- FindIntegrationAnchors(object.list = seu_list, reference = c(1:4), normalization.method = "SCT",
                                  anchor.features = features, dims = 1:50)
seu_int2 <- IntegrateData(anchorset = anchors, dims = 1:50)

rm(anchors, seu_list)

seu_int2 <- ScaleData(seu_int2, verbose = FALSE)
seu_int2 <- RunPCA(seu_int2, verbose = FALSE)

# saveRDS(seu_int, file = "RDSfiles/seu_int.RDS")
# seu_int <- readRDS(file = "RDSfiles/seu_int.RDS")

seu_int2 <- RunUMAP(seu_int2, dims = 1:50)
seu_int2 <- FindNeighbors(seu_int2, reduction = "pca", dims = 1:50)
seu_int2 <- FindClusters(seu_int2, resolution = 0.5)
DimPlot(seu_int2, label = TRUE, repel = TRUE) + NoLegend() + NoAxes()
DimPlot(seu_int2, group.by = "sample_ID") + NoAxes()

saveRDS(seu_int2, file = "RDSfiles/seu_int2.RDS")

# Check which cluster is epithelial, immune or stroma ----
seu_int2 <- readRDS(file = "RDSfiles/seu_int2.RDS")

DefaultAssay(seu_int2) <- "SCT"

FeaturePlot(seu_int2,features = "EPCAM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "MUC5AC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "TFF1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "LIPF", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "PGA3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "REG4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "TFF3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "CHGA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "TRPM5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "FN1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "RGS5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "NOTCH3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "LUM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "DCN", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "PLVAP", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "ACKR1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "PTPRC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "CD68", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "CD3D", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "CD79A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "CD8A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "IL2RA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "KLRD1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "GZMK", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "MS4A1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "TNFRSF17", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "KIT", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "PLD4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "CD163", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu_int2,features = "LRG1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "CD38", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "ENG", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "MUC6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "CD44", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "CLU", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int2,features = "S100A4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# subset epithelial, stromal and immune clusters
seu_epi2 <- subset(seu_int2, idents = c(3, 5, 9, 12, 15))
seu_st2 <- subset(seu_int2, idents = c(6, 7, 11, 16))
seu_imm2 <- subset(seu_int2, idents = c(0:2, 4, 8, 10, 13, 14, 17))

saveRDS(seu_epi2, file = "RDSfiles/seu_epi2.RDS")
saveRDS(seu_st2, file = "RDSfiles/seu_st2.RDS")
saveRDS(seu_imm2, file = "RDSfiles/seu_imm2.RDS")
