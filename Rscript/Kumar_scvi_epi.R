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
library(DESeq2)
library(fgsea)


# re-integrate the epithelial subset and then clean it---- 
seu_epi <- readRDS("RDSfiles/seu_epi_scvi.RDS")

seu_epi <-DietSeurat(seu_epi) 
seu_epi <- JoinLayers(seu_epi)
seu_epi[["RNA"]]$data <- NULL
seu_epi[["RNA"]]$scale.data <- NULL
seu_epi[["RNA"]] <- split(seu_epi[["RNA"]], f = seu_epi$Patient.ID)
seu_epi #sample 11 and 31 are missing (probably there are no epithelial cells in those samples)

seu_epi <- NormalizeData(seu_epi, verbose = FALSE)
seu_epi <- FindVariableFeatures(seu_epi, verbose = FALSE)
seu_epi <- ScaleData(seu_epi, verbose = FALSE)
seu_epi <- RunPCA(seu_epi, npcs = 30, verbose = FALSE)

seu_epi <- IntegrateLayers(
  object = seu_epi, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/usr/local/Caskroom/miniforge/base/envs/scvi-env2",
  verbose = TRUE
)

seu_epi <- FindNeighbors(seu_epi, reduction = "integrated.scvi", dims = 1:30, verbose = FALSE)
seu_epi <- FindClusters(seu_epi, resolution = 1, cluster.name = "scvi_clusters", verbose = FALSE)
seu_epi <- RunUMAP(seu_epi, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", verbose = FALSE)

DimPlot(seu_epi, label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu_epi, group.by = "sample_ID") + NoAxes()
DimPlot(seu_epi, group.by = "Tumor") + NoAxes()

FeaturePlot(seu_epi,features = "EPCAM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "PTPRC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "FN1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "MUC5AC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "MUC6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "LIPF", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "ATP4B", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "TFF3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "CHGA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "RGS5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "DCN", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "PLVAP", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "CD79A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "CD8A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "IL2RA", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "KLRD1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "MS4A1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "TNFRSF17", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "KIT", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu_epi,features = "LRG1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "CD38", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "ENG", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

saveRDS(seu_epi, file = "RDSfiles/seu_epi_scvi.RDS")

# remove clusters contaminated with stroma or immune
seu_epi <- subset(seu_epi, idents = c(0:16,18:25))
# repeat clustering and check feature plots until it is clean

saveRDS(seu_epi, file = "RDSfiles/seu_epi_scvi_cleaned.RDS")

# analysis only in tumor samples----
seu_epi_t <- subset(seu_epi, subset = Tumor == "Tumor")
Idents(seu_epi_t) <- "sample_ID"
seu_epi_t <- JoinLayers(seu_epi_t)

# rename idents according to differentiation markers----
Idents(seu_epi_t) <- "seurat_clusters"
DimPlot(seu_epi_t, label = TRUE, repel = TRUE) + NoLegend() + NoAxes()
Idents(seu_epi_t) <- "seurat_clusters"
seu_epi_t <- RenameIdents(seu_epi_t, 
                          `0` = "Intestinal", `1` = "Pit", `2` = "Intestinal", `3` = "Pit", `4` = "Chief", `5` = "Intestinal",`6` = "Chief", `7` = "Intestinal", `8` = "Pit", `9` = "Neck", 
                          `10` = "Pit", `11` = "Intestinal", `12` = "Pit", `13` = "Endocrine", `14` = "Pit", `15` = "Chief",`16` = "Intestinal", `17` = "Parietal", `18` = "Intestinal", `19` = "Endocrine", 
                          `20` = "Intestinal", `21` = "Pit", `22` = "Intestinal", `23` = "Endocrine", `24` = "Intestinal", `25` = "Endocrine",`26` = "Pit")
DimPlot(seu_epi_t, label = TRUE, repel = TRUE) + NoLegend() + NoAxes()
seu_epi_t$cell_type <- Idents(seu_epi_t)
seu_epi_t$cell_type <- factor(x = seu_epi_t$cell_type, levels = c("Pit", "Neck", "Chief", "Parietal", "Endocrine", "Intestinal"))
saveRDS(seu_epi_t, file = "RDSfiles/seu_epi_t.RDS")

# Rank samples by MUC6----
avg_epi_t <- AverageExpression(seu_epi_t, return.seurat = F, slot = "data", assays = "RNA", group.by = c("sample_ID")) # use "data" to rank (not "count")
write.table(as.matrix(avg_epi_t$RNA), "results/GeneExpression/avg_epi_t.txt", sep="\t",col.names = T, row.names = T)  # check it in excel, rank samples by MUC6 levels
Idents(seu_epi_t) <- "sample_ID"
table(seu_epi_t$sample_ID) # sample_ID 15 has only 13. remove that sample
seu_epi_t <- subset(seu_epi_t, idents = 15, invert = T)
MUC6_H <- WhichCells(seu_epi_t, ident = c(40,	24,	14,	16,	36,	29,	26,	5,	19,	38,	10,	8,	18,	2))
MUC6_L <- WhichCells(seu_epi_t, ident = c(22,	7,	27,	12,	33,	3,	13,	17,	34,	32,	30,	20,	28,	39))

seu_epi_t <- SetIdent(seu_epi_t, cells = MUC6_H, value = "MUC6_H")
seu_epi_t$MUC6_level <- Idents(seu_epi_t)
seu_epi_t <- SetIdent(seu_epi_t, cells = MUC6_L, value = "MUC6_L")
seu_epi_t$MUC6_level <- Idents(seu_epi_t)

saveRDS(seu_epi_t, file = "RDSfiles/seu_epi_t.RDS")

# seu_epi_t <- FindClusters(seu_epi_t, resolution = 0.1, cluster.name = "scvi_clusters", verbose = FALSE)
# DimPlot(seu_epi_t, label = F, repel = TRUE) + NoAxes()
DimPlot(seu_epi_t, group.by = "MUC6_level") + NoAxes()
DimPlot(seu_epi_t, group.by = "sample_ID") + NoAxes() + guides(color = guide_legend(override.aes = list(size=3), ncol=3))
# VlnPlot(seu_epi_t, features = c("MUC6", "GOLPH3", "TFE3", "EGFR"), group.by = "MUC6_level", pt.size = 0.0)

# pseudobulk analysis for DEseq2----
seu_epi_t <- JoinLayers(seu_epi_t)
bulk <- AggregateExpression(seu_epi_t, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("sample_ID", "MUC6_level"))
tail(Cells(bulk))
bulk$sample_ID <- sapply(strsplit(Cells(bulk), split = "_"), "[", 1)
bulk$MUC6_level <- sapply(strsplit(Cells(bulk), split = "-"), "[", 2)
bulk$MUC6_level <- factor(x = bulk$MUC6_level, levels = c("L", "H"))
Idents(bulk) <- "MUC6_level"
de_markers <- FindMarkers(bulk, ident.1 = "L", ident.2 = "H", slot = "counts", test.use = "DESeq2")
de_markers$gene <- rownames(de_markers)
write.table(as.matrix(de_markers),"results/DEG//DESeq2_Seurat_FindMarkers_MUC6LvsH.txt", sep ="\t", col.names = T,row.names = F)
ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,"")), colour = "red", size = 3)

# do DEseq2 from matrix make rank for fgsea----
cts <- as.matrix(bulk[["RNA"]]$counts)
conditions <- sapply(strsplit(Cells(bulk), split = "-"), "[", 2)
conditions <- factor(conditions, levels = c("L", "H"))
coldata <- cbind(colnames(cts), conditions) %>% as.data.frame()
coldata$conditions <- conditions
rownames(coldata) <- coldata$V1

dds <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata, design = ~conditions)
dds <- DESeq(dds)
resultsNames(dds)

vst = varianceStabilizingTransformation(dds)
pcaData <- plotPCA(vst, intgroup = "conditions", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=conditions)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()

res <- results(dds, contrast = c("conditions", "L", "H")) %>% data.frame()
# res <- lfcShrink(dds, contrast = c("conditions", "L", "H"), type="ashr") %>% data.frame()
res <- rownames_to_column(res, var = "gene_name")
write.table(as.matrix(res),"results/DEG//DESeq2_MUC6LvsH_1007.txt", sep ="\t", col.names = T,row.names = F)
ggplot(res, aes(log2FoldChange, -log10(pvalue))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(padj < 0.01, gene_name,"")), colour = "red", size = 3)

res2 <- res %>% 
  dplyr::select(gene_name, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene_name) %>% 
  summarize(stat=mean(stat))
ranks <- deframe(res2)

# MUC6 top and bottom 25%----
Idents(seu_epi_t) <- "sample_ID"
seu_epi_t2 <-subset(seu_epi_t, idents = c(38, 14, 40, 19, 36, 24, 10, 17, 13, 32, 30, 20, 28, 39))
MUC6_Hq <- WhichCells(seu_epi_t2, ident = c(38,14,40,19,36,24,10))
MUC6_Lq <- WhichCells(seu_epi_t2, ident = c(17,13,32,30,20,28,39))
seu_epi_t2 <- SetIdent(seu_epi_t2, cells = MUC6_Hq, value = "MUC6_Hq")
seu_epi_t2$MUC6_level2 <- Idents(seu_epi_t2)
seu_epi_t2 <- SetIdent(seu_epi_t2, cells = MUC6_Lq, value = "MUC6_Lq")
seu_epi_t2$MUC6_level2 <- Idents(seu_epi_t2)

bulk2 <- AggregateExpression(seu_epi_t2, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("sample_ID", "MUC6_level2"))
tail(Cells(bulk2))
bulk2$sample_ID <- sapply(strsplit(Cells(bulk2), split = "_"), "[", 1)
bulk2$MUC6_level2 <- sapply(strsplit(Cells(bulk2), split = "-"), "[", 2)
bulk2$MUC6_level2 <- factor(x = bulk2$MUC6_level2, levels = c("Lq", "Hq"))
Idents(bulk2) <- "MUC6_level2"
de_markers2 <- FindMarkers(bulk2, ident.1 = "Lq", ident.2 = "Hq", slot = "counts", test.use = "DESeq2")
de_markers2$gene <- rownames(de_markers2)

ggplot(de_markers2, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,"")), colour = "red", size = 3)
VlnPlot(bulk2, features = c("MUC6", "GOLPH3", "TFE3", "EGFR", "PCSK2", "FSTL3"), split.by = "MUC6_level2", cols = c("#377eb8", "#e41a1c"))

# fGSEA----
# prepare gene sets
collections <- list()
collections$BIOCARTA <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "BIOCARTA")
collections$CGP <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
collections$HALLMARKS <- msigdbr(species = "Homo sapiens", category = "H")
collections$C6 <- msigdbr(species = "Homo sapiens", category = "C6")
collections <- lapply(collections, function(x) {
  out <- split(x = x$gene_symbol, f = x$gs_name)
})
# run fgsea
fgseaRes <- fgsea(pathways = collections$C6, stats = ranks, eps=0.0, minSize = 10, maxSize = 500)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.25)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

plotEnrichment(collections$C6[["EGFR_UP.V1_UP"]], ranks) + labs(title="C6_EGFR_UP.V1_UP")
plotEnrichment(collections$C6[["MEK_UP.V1_UP"]], ranks) + labs(title="C6_MEK_UP.V1_UP")
plotEnrichment(collections$C6[["PDGF_ERK_DN.V1_DN"]], ranks) + labs(title="C6_PDGF_ERK_DN.V1_DN")
plotEnrichment(collections$HALLMARKS[["HALLMARK_KRAS_SIGNALING_DN"]], ranks) + labs(title="HALLMARK_KRAS_SIGNALING_DN")
plotEnrichment(collections$HALLMARKS[["HALLMARK_KRAS_SIGNALING_UP"]], ranks) + labs(title="HALLMARK_KRAS_SIGNALING_UP")

write.table(fgseaRes[,-8],"results/DEG/fgseaRes_MUC6KO_UP_DN_MUC6LvsH.txt", sep ="\t", col.names = T,row.names = F)

# use DE genes from MUC6KO mouse RNAseq
MUC6KO_UP_DN <- gmtPathways("gene_set/MUC6KO_UP_DN_hs.gmt")
fgseaRes <- fgsea(pathways = MUC6KO_UP_DN, stats = ranks, eps=0.0, minSize = 10, maxSize = 500)
plotEnrichment(MUC6KO_UP_DN[["MUC6KO_UP"]], ranks) + labs(title="MUC6KO_UP")
plotEnrichment(MUC6KO_UP_DN[["MUC6KO_DN"]], ranks) + labs(title="MUC6KO_DN")

# only intestinal tumor----
seu_int <- subset(seu_epi_t, subset = Laurens == "Intestinal")
avg_int<-AverageExpression(seu_int, group.by = "sample_ID" ) 
write.table(as.matrix(avg_int$RNA), "results/GeneExpression/avg_int.txt", sep="\t",col.names = T, row.names = T) # check MUC6 levels in excel
Idents(seu_int) <- "sample_ID"
seu_int <- subset(seu_int, ident = 22, invert = T)
MUC6_H <- WhichCells(seu_int, ident = c(14,7,26,29,18,16))
MUC6_L <- WhichCells(seu_int, ident = c(34,33,17,13,32,39))
seu_int <- SetIdent(seu_int, cells = MUC6_H, value = "MUC6_H")
seu_int$MUC6_level <- Idents(seu_int)
seu_int <- SetIdent(seu_int, cells = MUC6_L, value = "MUC6_L")
seu_int$MUC6_level <- Idents(seu_int)

# pseudobulk analysis of intestinal type
bulk <- AggregateExpression(seu_int, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("sample_ID", "MUC6_level"))
tail(Cells(bulk))
bulk$sample_ID <- sapply(strsplit(Cells(bulk), split = "_"), "[", 1)
bulk$MUC6_level <- sapply(strsplit(Cells(bulk), split = "-"), "[", 2)
bulk$MUC6_level <- factor(x = bulk$MUC6_level, levels = c("L", "H"))
Idents(bulk) <- "MUC6_level"
de_markers <- FindMarkers(bulk, ident.1 = "L", ident.2 = "H", slot = "counts", test.use = "DESeq2")
de_markers$gene <- rownames(de_markers)

write.table(as.matrix(de_markers),"./Results/DEG//DESea2_intestinal_Seurat_FindMarkers_MUC6LvsH.txt", sep ="\t", col.names = T,row.names = T)

ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,"")), colour = "red", size = 3)
VlnPlot(bulk, features = c("MUC6", "GOLPH3", "TFE3", "EGFR", "ATP4A", "DEFA5"), split.by = "MUC6_level", cols = c("#377eb8", "#e41a1c"))

# do DEseq2 from matrix make rank for fgsea
cts <- as.matrix(bulk[["RNA"]]$counts)
conditions <- sapply(strsplit(Cells(bulk), split = "-"), "[", 2)
conditions <- factor(conditions, levels = c("L", "H"))
coldata <- cbind(colnames(cts), conditions) %>% as.data.frame()
coldata$conditions <- conditions
rownames(coldata) <- coldata$V1

dds <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata, design = ~conditions)
dds <- DESeq(dds)

vst = varianceStabilizingTransformation(dds)
pcaData <- plotPCA(vst, intgroup = "conditions", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=conditions)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()

res <- results(dds, contrast = c("conditions", "L", "H")) %>% data.frame()
res <- rownames_to_column(res, var = "gene_name")
write.table(as.matrix(res),"results/DEG//DESeq2_Intestinal_MUC6LvsH.txt", sep ="\t", col.names = T,row.names = F)
ggplot(res, aes(log2FoldChange, -log10(pvalue))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(padj < 0.01, gene_name,"")), colour = "red", size = 3)

res2 <- res %>% 
  dplyr::select(gene_name, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene_name) %>% 
  summarize(stat=mean(stat))
ranks <- deframe(res2)

# run fgsea
fgseaRes <- fgsea(pathways = collections$C6, stats = ranks, eps=0.0, minSize = 10, maxSize = 500)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.25)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

plotEnrichment(collections$C6[["EGFR_UP.V1_UP"]], ranks) + labs(title="C6_EGFR_UP.V1_UP")
plotEnrichment(collections$C6[["MEK_UP.V1_UP"]], ranks) + labs(title="C6_MEK_UP.V1_UP")
plotEnrichment(collections$C6[["PDGF_ERK_DN.V1_DN"]], ranks) + labs(title="C6_PDGF_ERK_DN.V1_DN")
plotEnrichment(collections$HALLMARKS[["HALLMARK_KRAS_SIGNALING_DN"]], ranks) + labs(title="HALLMARK_KRAS_SIGNALING_DN")
plotEnrichment(collections$HALLMARKS[["HALLMARK_KRAS_SIGNALING_UP"]], ranks) + labs(title="HALLMARK_KRAS_SIGNALING_UP")

write.table(as.matrix(fgseaRes),"results/DEG//fgseaRes_intestinal_C6_MUC6LvsH.txt", sep ="\t", col.names = T,row.names = T)

# use DE genes from MUC6KO mouse RNAseq
MUC6KO_UP_DN <- gmtPathways("gene_set/MUC6KO_UP_DN_hs.gmt")
fgseaRes <- fgsea(pathways = MUC6KO_UP_DN, stats = ranks, eps=0.0, minSize = 10, maxSize = 500)
plotEnrichment(MUC6KO_UP_DN[["MUC6KO_UP"]], ranks) + labs(title="MUC6KO_UP")
plotEnrichment(MUC6KO_UP_DN[["MUC6KO_DN"]], ranks) + labs(title="MUC6KO_DN")
write.table(as.matrix(fgseaRes),"results/DEG//fgseaRes_intestinal_MUC6KO_UP_DN_MUC6LvsH.txt", sep ="\t", col.names = T,row.names = T)

# MKI67 positive cells----
seu_epi_t <- readRDS("RDSfiles/seu_epi_t.RDS")
Idents(seu_epi_t) <- "seurat_clusters"
FeaturePlot(seu_epi_t,features = "MKI67", cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
seu_epi_t1 <- subset(seu_epi_t, subset = MKI67 > 0.0)
DimPlot(seu_epi_t1)+ NoAxes()

# FeatureScatter
# seu_epi_t1 <- JoinLayers(seu_epi_t1)
# FeatureScatter(object = seu_epi_t1, feature1 = 'MUC6', feature2 = 'GOLPH3')
# FeatureScatter(object = seu_epi_t1, feature1 = 'MUC6', feature2 = 'TFE3')
# FeatureScatter(object = seu_epi_t1, feature1 = 'MUC6', feature2 = 'EGFR')

# calculate correlations with MUC6
# mtx<-seu_epi_t1[["RNA"]]$data
# mtx<-as.matrix(mtx)
# MUC6<-as.numeric(mtx["MUC6",])
# cors<-apply(mtx,1,function(x){cor(MUC6,x)})
# cors<-as.data.frame(cors)
# write.table(cors,"results/GeneExpression///MKI67pos_correlation_w_MUC6.txt", sep ="\t", col.names = T,row.names = T)

# re-integrate
seu_epi_t1 <-DietSeurat(seu_epi_t1) 
seu_epi_t1 <- JoinLayers(seu_epi_t1)
seu_epi_t1[["RNA"]]$data <- NULL
seu_epi_t1[["RNA"]]$scale.data <- NULL
seu_epi_t1[["RNA"]] <- split(seu_epi_t1[["RNA"]], f = seu_epi_t1$Patient.ID)
seu_epi_t1

seu_epi_t1 <- NormalizeData(seu_epi_t1, verbose = FALSE)
seu_epi_t1 <- FindVariableFeatures(seu_epi_t1, verbose = FALSE)
seu_epi_t1 <- ScaleData(seu_epi_t1, verbose = FALSE)
seu_epi_t1 <- RunPCA(seu_epi_t1, npcs = 30, verbose = FALSE)

seu_epi_t1 <- IntegrateLayers(
  object = seu_epi_t1, method = scVIIntegration,
  new.reduction = "integrated.scvi",
  conda_env = "/Users/kinoshitahiroto/miniconda3/envs/scvi-env",
  verbose = TRUE
)

seu_epi_t1 <- FindNeighbors(seu_epi_t1, reduction = "integrated.scvi", dims = 1:30, verbose = FALSE)
seu_epi_t1 <- FindClusters(seu_epi_t1, resolution = 1, cluster.name = "scvi_clusters", verbose = FALSE)
seu_epi_t1 <- RunUMAP(seu_epi_t1, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", verbose = FALSE)

DimPlot(seu_epi_t1, label = TRUE, repel = TRUE) + NoAxes()
DimPlot(seu_epi_t1, group.by = "sample_ID") + NoAxes()

FeaturePlot(seu_epi_t1,features = "MUC6", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi_t1,features = "GOLPH3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi_t1,features = "TFE3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi_t1,features = "EGFR", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi_t1,features = "MUC5AC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi_t1,features = "TFF1", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi_t1,features = "TFF2", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi_t1,features = "LIPF", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi_t1,features = "TFF3", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi_t1,features = "MUC2", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi_t1,features = "MUC4", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi_t1,features = "CDH17", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi_t1,features = "CD44", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi_t1,features = "CLU", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi_t1,features = "GKN3P", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

saveRDS(seu_epi_t1, "2023_Kumar_imacpro//seu_epi_t1.RDS")

VlnPlot(seu_epi_t1, features = c("MUC5AC", "MUC6", "TFF2", "TFF3", "MUC2", "CDH17", "CD44", "CLU", "GOLPH3", "MKI67"))

# pseudobulk
avg_epi_t1<-AverageExpression(seu_epi_t1, group.by = "sample_ID" ) 
write.table(as.matrix(avg_epi_t1$RNA), "2023_Kumar_imacpro/avg_epi_t1.txt", sep="\t",col.names = T, row.names = T)  # check it in excel
Idents(seu_epi_t) <- "sample_ID"
MUC6_H <- WhichCells(seu_epi_t, ident = c(2,40,14,8,27,38,36,5,29,7,22,26,17,24,10))
MUC6_L <- WhichCells(seu_epi_t, ident = c(16,13,34,19,33,12,30,18,32,28,3,15,20,39))

seu_epi_t1 <- SetIdent(seu_epi_t1, cells = MUC6_H, value = "MUC6_H")
seu_epi_t1$MUC6_level <- Idents(seu_epi_t1)
seu_epi_t1 <- SetIdent(seu_epi_t1, cells = MUC6_L, value = "MUC6_L")
seu_epi_t1$MUC6_level <- Idents(seu_epi_t1)
DimPlot(seu_epi_t1, group.by = "MUC6_level") + NoAxes()

bulk <- AggregateExpression(seu_epi_t1, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("sample_ID", "MUC6_level"))
tail(Cells(bulk))
bulk$sample_ID <- sapply(strsplit(Cells(bulk), split = "_"), "[", 1)
bulk$MUC6_level <- sapply(strsplit(Cells(bulk), split = "-"), "[", 2)
bulk$MUC6_level <- factor(x = bulk$MUC6_level, levels = c("L", "H"))
Idents(bulk) <- "MUC6_level"
de_markers <- FindMarkers(bulk, ident.1 = "L", ident.2 = "H", slot = "counts", test.use = "DESeq2")
de_markers$gene <- rownames(de_markers)

write.table(as.matrix(de_markers),"./Results/DEG//DESeq2_MKI67pos_MUC6LvsH.txt", sep ="\t", col.names = T,row.names = T)

ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,"")), colour = "red", size = 3)
VlnPlot(bulk, features = c("MUC6", "GOLPH3", "TFE3", "EGFR", "PCSK2", "FSTL3"), split.by = "MUC6_level", cols = c("#377eb8", "#e41a1c"))

# prepare gene sets
collections <- list()
collections$BIOCARTA <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "BIOCARTA")
collections$CGP <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
collections$HALLMARKS <- msigdbr(species = "Homo sapiens", category = "H")
collections$C6 <- msigdbr(species = "Homo sapiens", category = "C6")
collections <- lapply(collections, function(x) {
  out <- split(x = x$gene_symbol, f = x$gs_name)
})
# run fgsea
ranks <- de_markers$avg_log2FC
names(ranks) <- rownames(de_markers)
fgseaRes <- fgsea(pathways = collections$HALLMARKS, stats = ranks, eps=0.0, minSize = 10, maxSize = 500)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.25)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

# Rank samples by GOLPH3----
Idents(seu_epi_t) <- "sample_ID"
# check averaged gene expression in excel
GOLPH3_H <- WhichCells(seu_epi_t, ident = c(24,	39,	14,	18,	13,	2,	20,	8,	22,	7,	12,	19,	28,	36))
GOLPH3_L <- WhichCells(seu_epi_t, ident = c(33,	38,	40,	16,	17,	3,	34,	27,	5,	26,	10,	32,	30,	29))
seu_epi_t <- SetIdent(seu_epi_t, cells = GOLPH3_H, value = "GOLPH3_H")
seu_epi_t$GOLPH3_level <- Idents(seu_epi_t)
seu_epi_t <- SetIdent(seu_epi_t, cells = GOLPH3_L, value = "GOLPH3_L")
seu_epi_t$GOLPH3_level <- Idents(seu_epi_t)
seu_epi_t$GOLPH3_level <- factor(seu_epi_t$GOLPH3_level, levels = c("GOLPH3_H", "GOLPH3_L"))
DimPlot(seu_epi_t, group.by = "GOLPH3_level") + NoAxes()

# pseudobulk
bulk <- AggregateExpression(seu_epi_t, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("sample_ID", "GOLPH3_level"))
tail(Cells(bulk))
bulk$sample_ID <- sapply(strsplit(Cells(bulk), split = "_"), "[", 1)
bulk$GOLPH3_level <- sapply(strsplit(Cells(bulk), split = "-"), "[", 2)
bulk$GOLPH3_level <- factor(x = bulk$GOLPH3_level, levels = c("H", "L"))
Idents(bulk) <- "GOLPH3_level"
de_markers <- FindMarkers(bulk, ident.1 = "H", ident.2 = "L", slot = "counts", test.use = "DESeq2")
de_markers$gene <- rownames(de_markers)

write.table(as.matrix(de_markers),"./Results/DEG//DESeq2_Seurat_FindMarkers_GOLPH3LvsH.txt", sep ="\t", col.names = T,row.names = T)

ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,"")), colour = "red", size = 3)

# do DEseq2 from matrix make rank for fgsea
cts <- as.matrix(bulk[["RNA"]]$counts)
conditions <- sapply(strsplit(Cells(bulk), split = "-"), "[", 2)
conditions <- factor(conditions, levels = c("H", "L"))
coldata <- cbind(colnames(cts), conditions) %>% as.data.frame()
coldata$conditions <- conditions
rownames(coldata) <- coldata$V1

dds <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata, design = ~conditions)
dds <- DESeq(dds)

vst = varianceStabilizingTransformation(dds)
pcaData <- plotPCA(vst, intgroup = "conditions", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=conditions)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()

res <- results(dds, contrast = c("conditions", "H", "L")) %>% data.frame()
res <- rownames_to_column(res, var = "gene_name")
write.table(as.matrix(res),"results/DEG//DESeq2_GOLPH3LvsH.txt", sep ="\t", col.names = T,row.names = F)
ggplot(res, aes(log2FoldChange, -log10(pvalue))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(padj < 0.01, gene_name,"")), colour = "red", size = 3)

res2 <- res %>% 
  dplyr::select(gene_name, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene_name) %>% 
  summarize(stat=mean(stat))
ranks <- deframe(res2)

# run fgsea
fgseaRes <- fgsea(pathways = collections$C6, stats = ranks, eps=0.0, minSize = 10, maxSize = 500)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.25)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

plotEnrichment(collections$C6[["EGFR_UP.V1_UP"]], ranks) + labs(title="C6_EGFR_UP.V1_UP")
plotEnrichment(collections$C6[["MEK_UP.V1_UP"]], ranks) + labs(title="C6_MEK_UP.V1_UP")
plotEnrichment(collections$C6[["PDGF_ERK_DN.V1_DN"]], ranks) + labs(title="C6_PDGF_ERK_DN.V1_DN")
plotEnrichment(collections$HALLMARKS[["HALLMARK_KRAS_SIGNALING_DN"]], ranks) + labs(title="HALLMARK_KRAS_SIGNALING_DN")
plotEnrichment(collections$HALLMARKS[["HALLMARK_KRAS_SIGNALING_UP"]], ranks) + labs(title="HALLMARK_KRAS_SIGNALING_UP")

write.table(as.matrix(fgseaRes),"results/DEG//fgseaRes_C6_GOLPH3LvsH.txt", sep ="\t", col.names = T,row.names = T)

# use DE genes from MUC6KO mouse RNAseq
MUC6KO_UP_DN <- gmtPathways("gene_set/MUC6KO_UP_DN_hs.gmt")
fgseaRes <- fgsea(pathways = MUC6KO_UP_DN, stats = ranks, eps=0.0, minSize = 10, maxSize = 500)
plotEnrichment(MUC6KO_UP_DN[["MUC6KO_UP"]], ranks) + labs(title="MUC6KO_UP")
plotEnrichment(MUC6KO_UP_DN[["MUC6KO_DN"]], ranks) + labs(title="MUC6KO_DN")

write.table(as.matrix(fgseaRes),"results/DEG//fgseaRes_MUC6KO_UP_DN_GOLPH3LvsH.txt", sep ="\t", col.names = T,row.names = T)

# Rank samples by ARID1A----
seu_epi_t <- JoinLayers(seu_epi_t)
Idents(seu_epi_t) <- "sample_ID"
ARID1A_H <- WhichCells(seu_epi_t, ident = c(13,	30,	2,	29,	24,	18,	14,	20,	34,	19,	36,	16,	17,	22))
ARID1A_L <- WhichCells(seu_epi_t, ident = c(38,	33,	3,	8,	27,	28,	12,	5,	40,	39,	7,	32,	26,	10))
seu_epi_t <- SetIdent(seu_epi_t, cells = ARID1A_H, value = "ARID1A_H")
seu_epi_t$ARID1A_level <- Idents(seu_epi_t)
seu_epi_t <- SetIdent(seu_epi_t, cells = ARID1A_L, value = "ARID1A_L")
seu_epi_t$ARID1A_level <- Idents(seu_epi_t)
seu_epi_t$ARID1A_level <- factor(seu_epi_t$ARID1A_level, levels = c("ARID1A_L", "ARID1A_H"))
DimPlot(seu_epi_t, group.by = "ARID1A_level") + NoAxes()
saveRDS(seu_epi_t, "RDSfiles/seu_epi_t.RDS")

# pseudobulk following the Seurat vignette
bulk <- AggregateExpression(seu_epi_t, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("sample_ID", "ARID1A_level"))
tail(Cells(bulk))
bulk$sample_ID <- sapply(strsplit(Cells(bulk), split = "_"), "[", 1)
bulk$ARID1A_level <- sapply(strsplit(Cells(bulk), split = "-"), "[", 2)
bulk$ARID1A_level <- factor(x = bulk$ARID1A_level, levels = c("L", "H"))
Idents(bulk) <- "ARID1A_level"
# de_markers <- FindMarkers(bulk, ident.1 = "L", ident.2 = "H", slot = "counts", test.use = "DESeq2")
# de_markers$gene <- rownames(de_markers)
# write.table(as.matrix(de_markers),"./Results/DEG//DESeq2_Seurat_FindMarders_ARID1ALvsH.txt", sep ="\t", col.names = T,row.names = T)
# ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
#   ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,"")), colour = "red", size = 3)

# do DEseq2 from matrix make rank for fgsea
cts <- as.matrix(bulk[["RNA"]]$counts)
conditions <- sapply(strsplit(Cells(bulk), split = "-"), "[", 2)
conditions <- factor(conditions, levels = c("L", "H"))
coldata <- cbind(colnames(cts), conditions) %>% as.data.frame()
coldata$conditions <- conditions
rownames(coldata) <- coldata$V1

dds <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata, design = ~conditions)
dds <- DESeq(dds)

vst = varianceStabilizingTransformation(dds)
pcaData <- plotPCA(vst, intgroup = "conditions", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=conditions)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw()

res <- results(dds, contrast = c("conditions", "L", "H")) %>% data.frame()
res <- rownames_to_column(res, var = "gene_name")
write.table(as.matrix(res),"results/DEG//DESeq2_ARID1ALvsH.txt", sep ="\t", col.names = T,row.names = F)
ggplot(res, aes(log2FoldChange, -log10(pvalue))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(padj < 0.01, gene_name,"")), colour = "red", size = 3)

res2 <- res %>% 
  dplyr::select(gene_name, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene_name) %>% 
  summarize(stat=mean(stat))
ranks <- deframe(res2)

# run fgsea
fgseaRes <- fgsea(pathways = collections$C6, stats = ranks, eps=0.0, minSize = 10, maxSize = 500)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.25)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

write.table(fgseaRes[,-8],"results/DEG/fgseaRes_C6_ARID1ALvsH.txt", sep ="\t", col.names = T,row.names = F)

plotEnrichment(collections$C6[["AKT_UP.V1_DN"]], ranks) + labs(title="C6_AKT_UP.V1_DN")
plotEnrichment(collections$C6[["AKT_UP.V1_UP"]], ranks) + labs(title="C6_AKT_UP.V1_UP")

# fgsea with several selected gene sets
# use DE genes from ARID1AKO mouse RNAseq
ARID1AKO_UP_DN <- gmtPathways("gene_set/ARID1AKO_UP_DN_hs.gmt")
fgseaRes <- fgsea(pathways = ARID1AKO_UP_DN, stats = ranks, eps=0.0, minSize = 10, maxSize = 500)
plotEnrichment(ARID1AKO_UP_DN[["ARID1AKO_UP"]], ranks) + labs(title="ARID1AKO_UP")
plotEnrichment(ARID1AKO_UP_DN[["ARID1AKO_DN"]], ranks) + labs(title="ARID1AKO_DN")
write.table(fgseaRes[,-8],"results/DEG/fgseaRes_ARID1AKO_UP_DN_ARID1ALvsH.txt", sep ="\t", col.names = T,row.names = F)

KEGG_ASTHMA <- gmtPathways("gene_set/KEGG_ASTHMA.v2023.2.Hs.gmt")
fgseaRes <- fgsea(pathways = KEGG_ASTHMA, stats = ranks, eps=0.0, minSize = 10, maxSize = 500)
plotEnrichment(KEGG_ASTHMA[["KEGG_ASTHMA"]], ranks) + labs(title="KEGG_ASTHMA")
write.table(fgseaRes[,-8],"results/DEG/fgseaRes_KEGG_ASTHMA_ARID1ALvsH.txt", sep ="\t", col.names = T,row.names = F)


# bulk of data slot for vln plot
bulk2 <- AggregateExpression(seu_epi_t, return.seurat = T, slot = "data", assays = "RNA", group.by = c("sample_ID", "ARID1A_level"))
tail(Cells(bulk2))
bulk2$sample_ID <- sapply(strsplit(Cells(bulk2), split = "_"), "[", 1)
bulk2$ARID1A_level <- sapply(strsplit(Cells(bulk2), split = "-"), "[", 2)
bulk2$ARID1A_level <- factor(x = bulk2$ARID1A_level, levels = c("L", "H"))
Idents(bulk) <- "ARID1A_level"
VlnPlot(bulk2, features = c("IL33", "BTC"), group.by = "ARID1A_level")

FeaturePlot(seu_epi_t,features = "IL33", cols = c("lightgrey","darkred"), label = F, repel = TRUE) + NoAxes() + NoLegend()
