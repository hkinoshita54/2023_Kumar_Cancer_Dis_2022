library(tidyverse)
library(data.table)
library(Matrix)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(Seurat)
options(Seurat.object.assay.version = "v5")
library(SeuratWrappers)
library(msigdbr)
library(fgsea)
library(speckle)
library(miloR)
library(SingleCellExperiment)
library(scater)

# DA analysis by miloR----
seu_imm_t <- JoinLayers(seu_imm_t)
options(Seurat.object.assay.version = "v3")
seu_imm_t[["RNA"]] <- as(object = seu_imm_t[["RNA"]], Class = "Assay")
sce_imm_t <- as.SingleCellExperiment(seu_imm_t)
milo_imm_t <- Milo(sce_imm_t)
milo_imm_t <- buildGraph(milo_imm_t, k = 100, d = 30, reduced.dim = "PCA")
milo_imm_t <- makeNhoods(milo_imm_t, prop = 0.2, k = 100, d=30, refined = TRUE, reduced_dims = "PCA")
plotNhoodSizeHist(milo_imm_t)

milo_imm_t <- countCells(milo_imm_t, meta.data = data.frame(colData(milo_imm_t)), samples="sample_ID")
head(nhoodCounts(milo_imm_t))

design <- data.frame(colData(milo_imm_t))[,c("sample_ID", "IL33_level")]
design <- distinct(design)
rownames(design) <- design$sample_ID
design <- design[colnames(nhoodCounts(milo_imm_t)), , drop=FALSE]
design

milo_imm_t <- calcNhoodDistance(milo_imm_t, d=30, reduced.dim = "PCA")
rownames(design) <- design$sample_ID
contrast.1 <- c("IL33_levelIL33_H - IL33_levelIL33_L")
da_results <- testNhoods(milo_imm_t, design = ~ 0 + IL33_level, design.df = design, model.contrasts = contrast.1)
da_results %>%
  arrange(- SpatialFDR) %>%
  head() 

milo_imm_t <- buildNhoodGraph(milo_imm_t)

ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) ## Mark significance threshold (10% FDR)


milo_imm_t <- buildNhoodGraph(milo_imm_t)

## Plot single-cell UMAP
umap_pl <- plotReducedDim(milo_imm_t, dimred = "UMAP.SCVI", colour_by="IL33_level", text_by = "cellSubtype", 
                          text_size = 3, point_size=0.5) +
  guides(fill="none")

## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(milo_imm_t, da_results, layout="UMAP.SCVI",alpha=0.1) 

umap_pl + nh_graph_pl + plot_layout(guides="collect")

da_results <- annotateNhoods(milo_imm_t, da_results, coldata_col = "cellSubtype")
head(da_results)
write.table(da_results, file = "results/DA//miloR_res_imm_t_DF_annotated_IL33HvsL.txt", sep ="\t", col.names = T,row.names = F)
plotDAbeeswarm(da_results, group.by = "cellSubtype")

# DA analysis by speckle----
propres <- propeller(clusters = seu_imm_t$cellSubtype, sample = seu_imm_t$sample_ID, group = seu_imm_t$IL33_level)
write.table(propres, file = "results/DA//propeller_imm_t_DF_annotated_IL33HvsL.txt", sep ="\t", col.names = T,row.names = F)
plotCellTypeProps(clusters = seu_imm_t$cellSubtype, sample = seu_imm_t$sample_ID)
