library(tidyverse)
library(data.table)
library(Matrix)
library(cowplot)
library(Seurat)
library(ggplot2)

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

# "integrating large datasets" pipline from Seurat vignettes----
rm(meta, meta_list, mtx, mtx_list, seu, files_list, i)
gc()
seu_list <- lapply(X = seu_list, FUN = function(x){
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = seu_list)
seu_list <- lapply(X = seu_list, FUN = function(x){
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verboose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = seu_list, reference = c(1:4), reduction = "rpca",
                                  dims = 1:50)
seu_int_2 <- IntegrateData(anchorset = anchors, dims = 1:50)

rm(anchors, seu_list)

seu_int_2 <- ScaleData(seu_int_2, verbose = FALSE)
seu_int_2 <- RunPCA(seu_int_2, verbose = FALSE)

# saveRDS(seu_int_2, file = "RDSfiles/seu_int_2.RDS")
# seu_int_2 <- readRDS(file = "RDSfiles/seu_int_2.RDS")

seu_int_2 <- RunUMAP(seu_int_2, dims = 1:30)
seu_int_2 <- FindNeighbors(seu_int_2, reduction = "pca", dims = 1:30)
seu_int_2 <- FindClusters(seu_int_2, resolution = 0.5)
DimPlot(seu_int_2, label = TRUE, repel = TRUE) + NoLegend()
DimPlot(seu_int_2, group.by = "sample_ID")

saveRDS(seu_int_2, file = "RDSfiles/seu_int_2.RDS")

# Check which cluster is epithelial, immune or stroma ----
seu_int_2 <- readRDS(file = "RDSfiles/seu_int_2.RDS")

DefaultAssay(seu_int_2) <- "RNA"
seu_int_2 <- NormalizeData(seu_int_2) %>% FindVariableFeatures() %>% ScaleData()
markers <- c("EPCAM", "CDH1", "MUC5AC", "TFF1", "LIPF", "PGA3", "REG4", "TFF3", "CHGA", "TRPM5", "FN1", "RGS5", "NOTCH3", "LUM", "DCN",
             "PLVAP", "ACKR1", "PTPRC", "CD68", "CD3D", "CD79A", "CD8A", "IL2RA", "KLRD1", "MS4A1", "TNFRSF17", "KIT", "PLD4", "CD163")

pdf(file = "plots/2_int_feature_plots.pdf")
# par(mfrow = c(5, 5))

for (i in seq_along(markers)) {
  plot <- FeaturePlot(seu_int_2,features = markers[i], cols = c("lightgrey","darkred"), 
                      label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
  print(plot)
}

dev.off()

FeaturePlot(seu_int_2,features = "PTGS1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# subset epithelial, stromal and immune clusters
seu_epi_2 <- subset(seu_int_2, idents = c(2, 7, 8, 9, 19, 21))
seu_st_2 <- subset(seu_int_2, idents = c(4, 5, 15, 18))
seu_imm_2 <- subset(seu_int_2, idents = c(0, 1, 3, 6, 10, 11, 12, 13, 14, 16, 17, 20))

saveRDS(seu_epi_2, file = "RDSfiles/seu_epi_2.RDS")
saveRDS(seu_st_2, file = "RDSfiles/seu_st_2.RDS")
saveRDS(seu_imm_2, file = "RDSfiles/seu_imm_2.RDS")
