library(tidyverse)
library(data.table)
library(cowplot)
library(Seurat)
library(ggplot2)

# make a mata_data table for each saample----
## patients_data.txt was made with excel from suppl. table1 and series_matrix
patients_data <- read.delim("meta_data/patients_data.txt", header = TRUE)
mtx_list <- list.files(path = "raw_data/GSE183904_RAW", pattern = "*.csv.gz", full.names = T) %>% 
  lapply(fread, header = TRUE) %>% 
  lapply(data.frame, row.names = 1)

meta_list <- list()
for (i in seq_along(mtx_list)) {
  meta <- data.frame(cell_ID = colnames(mtx_list[[i]]), sample_ID = i)
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

# simply merge and do SCTransform assuming the data are already integrated----
## requires large memory
# seu_merge <- merge(seu_list[[1]], seu_list[2:length(seu_list)], add.cell.ids = seq_along(seu_list))
# 
# seu_merge[["percent.mt"]] <- PercentageFeatureSet(seu_merge, pattern = "^MT-")
# VlnPlot(seu_merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
# 
# rm(meta, meta_list, mtx_list, patients_data, seu, seu_list)
# gc()
# 
# seu_merge <- SCTransform(seu_merge, method = "glmGamPoi", verbose = TRUE)
# seu_merge <- RunPCA(seu_merge, verbose = T)
# ElbowPlot(seu_merge)
# 
# saveRDS(seu_merge, file = "RDSfiles/seu_merge.RDS")
# 
# seu_merge <- RunUMAP(seu_merge, dims = 1:15, verbose = T)
# seu_merge <- FindNeighbors(seu_merge, dims = 1:15, verbose = T)
# seu_merge <- FindClusters(seu_merge, verbose = T, resolution = 0.5)
# 
# DimPlot(seu_merge, label = TRUE, repel = TRUE) + NoLegend()
# DimPlot(seu_merge, label = TRUE, repel = TRUE, group.by = "sample_ID")   

# "integrating large datasets" pipline from Seurat vignettes----
rm(meta, meta_list, mtx_list, seu)
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

anchors <- FindIntegrationAnchors(object.list = seu_list, reference = c(4:7), reduction = "rpca",
                                  dims = 1:50)
seu_int <- IntegrateData(anchorset = anchors, dims = 1:50)

rm(anchors, seu_list)

seu_int <- ScaleData(seu_int, verbose = FALSE)
seu_int <- RunPCA(seu_int, verbose = FALSE)

# saveRDS(seu_int, file = "RDSfiles/seu_int.RDS")
# seu_int <- readRDS(file = "RDSfiles/seu_int.RDS")

seu_int <- RunUMAP(seu_int, dims = 1:30)
seu_int <- FindNeighbors(seu_int, reduction = "pca", dims = 1:30)
seu_int <- FindClusters(seu_int, resolution = 0.5)
DimPlot(seu_int, label = TRUE, repel = TRUE) + NoLegend()
DimPlot(seu_int, group.by = "sample_ID")

saveRDS(seu_int, file = "RDSfiles/seu_int.RDS")

# Check which cluster is epithelial, immune or stroma ----
seu_int <- readRDS(file = "RDSfiles/seu_int.RDS")

DefaultAssay(seu_int) <- "RNA"
seu_int <- NormalizeData(seu_int) %>% FindVariableFeatures() %>% ScaleData()

markers <- c("EPCAM", "CDH1", "MUC5AC", "TFF1", "LIPF", "PGA3", "REG4", "TFF3", "CHGA", "TRPM5", "FN1", "RGS5", "NOTCH3", "LUM", "DCN",
             "PLVAP", "ACKR1", "PTPRC", "CD68", "CD3D", "CD79A", "CD8A", "IL2RA", "KLRD1", "MS4A1", "TNFRSF17", "KIT", "PLD4", "CD163")

pdf(file = "plots/1_int_feature_plots.pdf")
# par(mfrow = c(5, 5))

for (i in seq_along(markers)) {
  plot <- FeaturePlot(seu_int,features = markers[i], cols = c("lightgrey","darkred"), 
                      label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
  print(plot)
}

dev.off()


FeaturePlot(seu_int,features = "TMPRSS", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# subset epithelial, stromal and immune clusters----
seu_epi <- subset(seu_int, idents = c(3, 6, 11, 17))
seu_st <- subset(seu_int, idents = c(5, 7, 9, 15))
seu_imm <- subset(seu_int, idents = c(0:2, 4, 8, 10, 12:14, 16, 18))

saveRDS(seu_epi, file = "RDSfiles/seu_epi.RDS")
saveRDS(seu_st, file = "RDSfiles/seu_st.RDS")
saveRDS(seu_imm, file = "RDSfiles/seu_imm.RDS")


