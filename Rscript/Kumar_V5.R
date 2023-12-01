library(tidyverse)
library(data.table)
library(Matrix)
library(cowplot)
library(ggplot2)
library(patchwork)
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratWrappers)
options(future.globals.maxSize = 1e10)
library(harmony)

# Sys.setenv(RETICULATE_MINICONDA_PATH = '/Users/kinoshitahiroto/miniconda3/')
library(reticulate)
# use_miniconda('/Users/kinoshitahiroto/miniconda3/envs/scvi-env')

# make a mata_data table for each saample----
## patients_data.txt was made with excel from suppl. table1 and series_matrix
patients_data <- read.delim("meta_data/patients_data.txt", header = TRUE, sep = " ")
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
  meta <- data.frame(cell_ID = colnames(mtx_list[[i]]), sample_ID = i + 5)
  meta <- left_join(meta, patients_data)
  meta <- column_to_rownames(meta, var = "cell_ID")
  meta_list[[i]] <- meta
}

# Create Seurat object----
names(mtx_list) <- paste0("sample", seq_along(mtx_list))
meta_merge <- Reduce(rbind, meta_list)
seu <- CreateSeuratObject(counts = mtx_list, meta.data = meta_merge)

# "integrative analysis in Seurat V5" pipline from the vignettes----
rm(meta, meta_list, meta_merge, mtx, mtx_list, files_list, i)
gc()

seu$percent.mt <- PercentageFeatureSet(seu, pattern = "^MT-")
# VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
# seu <- subset(seu, subset = nCount_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 12)

# check cell number in each sample
# table(seu$sample_ID)
# ncell <- table(seu$sample_ID) %>% as.data.frame()
# ncell <- ncell[,2]
# patients_data$ncell <- ncell
# write_delim(patients_data, "meta_data/patients_data.txt")

# when useing SCTransform run below
# seu <- SCTransform(seu, method = "glmGamPoi", vars.to.regress = "percent.mt",
#                    vst.flavor = "v2", verbose = FALSE)

# the standard (lognormalization) pipeline from the vignette
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 30, verbose = FALSE)

saveRDS(seu, file = "RDSfiles/seu_V5_lognorm.RDS")

# Dimplot without integration
# seu <- FindNeighbors(seu, reduction = "pca", dims = 1:30, verbose = FALSE)
# seu <- FindClusters(seu, resolution = 0.8, cluster.name = "unintegrated_clusters", verbose = FALSE)
# seu <- RunUMAP(seu, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated", verbose = FALSE)
# DimPlot(seu, reduction = "umap.unintegrated", group.by = c("sample_ID", "unintegrated_clusters"))

seu <- readRDS(file = "RDSfiles/seu_V5_lognorm.RDS")

# Dimplot with integration
seu <- IntegrateLayers(
  object = seu, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  # reference = c(6, 7),
  conda_env = "/Users/kinoshitahiroto/miniconda3/envs/scvi-env",
  verbose = TRUE
)

# seu <- IntegrateLayers(
#   object = seu, method = HarmonyIntegration,
#   orig.reduction = "pca", new.reduction = "integrated.harmony",
#   normalization.method = "SCT",
#   # reference = c(6, 7),
#   verbose = FALSE
# )

seu <- FindNeighbors(seu, reduction = "integrated.scvi", dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 1, cluster.name = "scvi_clusters", verbose = FALSE)
seu <- RunUMAP(seu, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", verbose = FALSE)
# DimPlot(seu, reduction = "umap.scvi",
#         group.by = c("sample_ID", 
#                      # "unintegrated_clusters", 
#                      "scvi_clusters"),
#         combine = FALSE
# )
DimPlot(seu, label = TRUE, repel = TRUE) + NoLegend() + NoAxes()
DimPlot(seu, group.by = "sample_ID") + NoAxes()

saveRDS(seu, file = "RDSfiles/seu_V5_lognorm_scvi.RDS")

# Check which cluster is epithelial, immune or stroma ----
# seu <- readRDS(file = "RDSfiles/seu_V5_lognorm_rpca.RDS")

# DefaultAssay(seu) <- "RNA"
# seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()

FeaturePlot(seu,features = "EPCAM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "PTPRC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "FN1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "MUC5AC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "TFF1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "LIPF", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "PGA3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "REG4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "TFF3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "CHGA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "RGS5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "NOTCH3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "LUM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "DCN", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "PLVAP", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "ACKR1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "CD68", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "CD79A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "CD8A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "IL2RA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "KLRD1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "MS4A1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "TNFRSF17", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "KIT", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "PLD4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu,features = "CD163", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()


# subset epithelial, stromal and immune clusters
seu_epi <- subset(seu, idents = c(9, 13, 14, 22, 24, 27, 28, 32))
seu_st <- subset(seu, idents = c(6, 7, 10, 12, 15, 17, 20, 21, 25, 29, 30, 34))
seu_imm <- subset(seu, idents = c(0:5, 8, 11, 16, 18, 19, 23, 26, 31, 33))

saveRDS(seu_epi, file = "RDSfiles/seu_epi_v5.RDS")
saveRDS(seu_st, file = "RDSfiles/seu_st_v5.RDS")
saveRDS(seu_imm, file = "RDSfiles/seu_imm_v5.RDS")
