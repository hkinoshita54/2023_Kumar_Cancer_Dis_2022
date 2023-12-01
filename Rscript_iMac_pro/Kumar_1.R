library(tidyverse)
library(data.table)
library(cowplot)
library(Seurat)
library(ggplot2)
library(Matrix)
library(glmGamPoi)

# make a mata_data table for each saample----
## patients_data.txt was made with excel from suppl. table1 and series_matrix
patients_data <- read.delim("meta_data/patients_data.txt", header = TRUE)
files_list <- list.files(path = "raw_data/GSE183904_RAW", pattern = "*.csv.gz", full.names = T)
mtx_list <- list()
for (i in seq_along(files_list)) {
  mtx <- fread(files_list[[i]], header = TRUE)
  mtx <- data.frame(mtx, row.names = 1)
  mtx <- Matrix(as.matrix(mtx), sparse = TRUE)
  mtx_list[[i]] <- mtx
}

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
  x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = c("percent.mt"), verbose = TRUE)
})

features <- SelectIntegrationFeatures(object.list = seu_list, nfeatures = 3000)
seu_list <- PrepSCTIntegration(object.list = seu_list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = seu_list, reference = c(4:7), normalization.method = "SCT",
                                  anchor.features = features, dims = 1:50)
seu_int <- IntegrateData(anchorset = anchors, dims = 1:50)

rm(anchors, seu_list)

seu_int <- ScaleData(seu_int, verbose = FALSE)
seu_int <- RunPCA(seu_int, verbose = FALSE)

# saveRDS(seu_int, file = "RDSfiles/seu_int.RDS")
# seu_int <- readRDS(file = "RDSfiles/seu_int.RDS")

seu_int <- RunUMAP(seu_int, reduction = "pca", dims = 1:50)
seu_int <- FindNeighbors(seu_int, reduction = "pca", dims = 1:50)
seu_int <- FindClusters(seu_int, resolution = 0.5)
DimPlot(seu_int, label = TRUE, repel = TRUE) + NoLegend() + NoAxes()
DimPlot(seu_int, group.by = "sample_ID") + NoAxes()

saveRDS(seu_int, file = "RDSfiles/seu_int.RDS")

# Check which cluster is epithelial, immune or stroma ----
seu_int <- readRDS(file = "RDSfiles/seu_int.RDS")

DefaultAssay(seu_int) <- "SCT"

DefaultAssay(seu_int) <- "RNA"
seu_int <- NormalizeData(seu_int) %>% FindVariableFeatures() %>% ScaleData()

# markers <- c("EPCAM", "CDH1", "MUC5AC", "TFF1", "LIPF", "PGA3", "REG4", "TFF3", "CHGA", "TRPM5", "FN1", "RGS5", "NOTCH3", "LUM", "DCN",
#              "PLVAP", "ACKR1", "PTPRC", "CD68", "CD3D", "CD79A", "CD8A", "IL2RA", "KLRD1", "MS4A1", "TNFRSF17", "KIT", "PLD4", "CD163")

# pdf(file = "plots/1_int_feature_plots.pdf")
# par(mfrow = c(5, 5))
# 
# for (i in seq_along(markers)) {
#   plot <- FeaturePlot(seu_int,features = markers[i], cols = c("lightgrey","darkred"), 
#                       label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
#   print(plot)
# }
# 
# dev.off()

FeaturePlot(seu_int,features = "EPCAM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "MUC5AC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "TFF1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "LIPF", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "PGA3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "REG4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "TFF3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "CHGA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "TRPM5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "FN1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "RGS5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "NOTCH3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "LUM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "DCN", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "PLVAP", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "ACKR1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "PTPRC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "CD68", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "CD3D", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "CD79A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "CD8A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "IL2RA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "KLRD1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "GZMK", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "MS4A1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "TNFRSF17", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "KIT", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "PLD4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "CD163", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu_int,features = "LRG1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "CD38", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "ENG", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "MUC6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "CD44", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_int,features = "CLU", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()



# subset epithelial, stromal and immune clusters----
seu_epi <- subset(seu_int, idents = c(4, 5, 16, 18))
seu_st <- subset(seu_int, idents = c(6, 7, 8, 14))
seu_imm <- subset(seu_int, idents = c(0:3, 4, 9:13, 17, 19, 20))

saveRDS(seu_epi, file = "RDSfiles/seu_epi.RDS")
saveRDS(seu_st, file = "RDSfiles/seu_st.RDS")
saveRDS(seu_imm, file = "RDSfiles/seu_imm.RDS")


