library(tidyverse)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratWrappers)
library(reticulate)
options(future.globals.maxSize = 1e10)
library(msigdbr)
library(fgsea)
library(slingshot)


# sligshot
# https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/slingshot/slingshot.html

seu_epi_t <- readRDS("2023_Kumar_imacpro/seu_epi_t.RDS")
seu_epi_t <- JoinLayers(seu_epi_t)
dimred <- seu_epi_t@reductions$umap@cell.embeddings
clustering <- seu_epi_t$scvi_clusters
counts <- as.matrix(seu_epi_t[["RNA"]]$counts[rownames(seu_epi_t[["RNA"]]$scale.data), ])

set.seed(1)
lineages <- getLineages(data = dimred, clusterLabels = clustering)
lineages

pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))
par(mfrow = c(1, 2))
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
for (i in levels(clustering)) {
  text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
}
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
lines(lineages, lwd = 3, col = "black")
