library(tidyverse)
library(data.table)
library(Matrix)
library(cowplot)
library(ggplot2)
library(patchwork)
library(Seurat)
options(Seurat.object.assay.version = "v3")
library(SeuratWrappers)
library(DoubletFinder)
library(reticulate)
options(future.globals.maxSize = 1e10)

# make a mata_data table for each saample----
patients_data <- read.delim("meta_data/patients_data.txt", header = TRUE, sep = " ") # patients_data.txt was made with excel from suppl. table1 and series_matrix
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
  seu <- CreateSeuratObject(mtx_list[[i]], min.cells = 3, min.features = 200, meta.data = meta_list[[i]])
  seu_list[[i]] <- seu
}

# estimation of doublet rate is from 10x genomics web page: 0.8% for 1000 cells
# https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
ncell_recov <- lapply(mtx_list, ncol) %>% unlist()
doublet_rate <- 0.008/1000 * ncell_recov

# filter by percent.mt < 12----
for (i in seq_along(seu_list)) {
  seu_list[[i]]$percent.mt <- PercentageFeatureSet(seu_list[[i]], pattern = "^MT-")
  seu_list[[i]] <- subset(seu_list[[i]], subset = percent.mt < 12)
}

# VlnPlot(seu_list[[7]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 4, pt.size = 0)

rm(meta, meta_list, mtx, mtx_list, seu, files_list, i)

# loop through samples to find doublets----
for (i in 1:length(seu_list)) {
  # print the sample we are on
  print(paste0("Sample ",i))
  
  # Pre-process seurat object with standard seurat workflow
  seu <- NormalizeData(seu_list[[i]])
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu)
  seu <- RunUMAP(seu, dims = 1:10)
  seu <- FindNeighbors(object = seu, dims = 1:10)              
  seu <- FindClusters(object = seu, resolution = 0.1)
  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep_v3(seu, PCs = 1:10)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  annotations <- seu@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(doublet_rate[i] * nrow(seu@meta.data)) ## 
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  seu <- doubletFinder_v3(seu = seu, 
                          PCs = 1:10, 
                          pK = optimal.pk,
                          nExp = nExp.poi.adj)
  metadata <- seu@meta.data
  colnames(metadata)[ncol(metadata)] <- "doublet_finder"
  seu@meta.data <- metadata 
  
  # subset and save
  singlets <- subset(seu, doublet_finder == "Singlet")
  seu_list[[i]] <- singlets
  remove(singlets)
}

# merge seurat objects and convert to V5 object----
seu <- merge(seu_list[[1]], seu_list[2:length(seu_list)])
seu[["RNA"]] <- as(object = seu[["RNA"]], Class = "Assay5")
seu@meta.data <- seu@meta.data[1:22]

saveRDS(seu, file = "RDSfiles/seu_DF.RDS")
