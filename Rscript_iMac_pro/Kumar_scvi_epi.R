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
library(DESeq2)

# re-integrate epithelial subset and clean it to remove immune or stroma contaminatino----
# latest analysis if from Kumar_DF, in which mt<12 cells were selected
seu_epi <- readRDS("RDSfiles/seu_epi_scvi.RDS")

seu_epi <-DietSeurat(seu_epi) 
seu_epi <- JoinLayers(seu_epi)
seu_epi[["RNA"]]$data <- NULL
seu_epi[["RNA"]]$scale.data <- NULL
seu_epi[["RNA"]] <- split(seu_epi[["RNA"]], f = seu_epi$sample_ID)
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
seu_epi <- FindClusters(seu_epi, resolution = 1.2, cluster.name = "scvi_clusters", verbose = FALSE)
seu_epi <- RunUMAP(seu_epi, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi", verbose = FALSE)

DimPlot(seu_epi, label = T, repel = TRUE) + NoAxes()
DimPlot(seu_epi, group.by = "sample_ID") + NoAxes()
DimPlot(seu_epi, group.by = "Tumor") + NoAxes()

FeaturePlot(seu_epi,features = "EPCAM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "PTPRC", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "VIM", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

FeaturePlot(seu_epi,features = "MUC5AC", cols = c("lightgrey","darkred"), label = T, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "MUC6", cols = c("lightgrey","darkred"), label = T, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "LIPF", cols = c("lightgrey","darkred"), label = T, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "ATP4B", cols = c("lightgrey","darkred"), label = T, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "TFF3", cols = c("lightgrey","darkred"), label = T, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "CHGA", cols = c("lightgrey","darkred"), label = T, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "RGS5", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "DCN", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "PLVAP", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "CD79A", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "CD3D", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "CD68", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
FeaturePlot(seu_epi,features = "TNFRSF17", cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()

# remove clusters contaminated with stroma or immune
seu_epi <- subset(seu_epi, idents = c(12,14), invert = TRUE)
# repeat clustering and check feature plots until it is clean

saveRDS(seu_epi, file = "RDSfiles/seu_epi_DF.RDS")

# analysis only in tumor samples----
seu_epi_t <- subset(seu_epi, subset = Tumor == "Tumor")
table(seu_epi_t$sample_ID) #sample 3, 13 and 15 have <100 cells
Idents(seu_epi_t) <- "sample_ID"
seu_epi_t <- subset(seu_epi_t, idents = c(3,13,15), invert = TRUE)
seu_epi_t <- JoinLayers(seu_epi_t)

# rank samples by MUC6 level----
avg_epi_t<-AverageExpression(seu_epi_t, group.by = "sample_ID" ) 
write.table(as.matrix(avg_epi_t$RNA), "results/GeneExpression/avg_epi_t.txt", sep="\t",col.names = T, row.names = T)  # check it in excel
Idents(seu_epi_t) <- "sample_ID"
MUC6_H <- WhichCells(seu_epi_t, ident = c(38,	14,	40,	19,	24,	36,	5,	26,	10,	7,	18,	2,	16))
MUC6_L <- WhichCells(seu_epi_t, ident = c(22,	8,	27,	34,	12,	33,	29,	17,	30,	20,	32,	28,	39))

seu_epi_t <- SetIdent(seu_epi_t, cells = MUC6_H, value = "MUC6_H")
seu_epi_t$MUC6_level <- Idents(seu_epi_t)
seu_epi_t <- SetIdent(seu_epi_t, cells = MUC6_L, value = "MUC6_L")
seu_epi_t$MUC6_level <- Idents(seu_epi_t)
DimPlot(seu_epi_t, group.by = "MUC6_level") + NoAxes()

# pseudobulk analysis for DESeq2 following Seurat vignette----
bulk <- AggregateExpression(seu_epi_t, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("sample_ID", "MUC6_level"))
tail(Cells(bulk))
bulk$sample_ID <- sapply(strsplit(Cells(bulk), split = "_"), "[", 1)
bulk$MUC6_level <- sapply(strsplit(Cells(bulk), split = "-"), "[", 2)
bulk$MUC6_level <- factor(x = bulk$MUC6_level, levels = c("L", "H"))
Idents(bulk) <- "MUC6_level"
de_markers <- FindMarkers(bulk, ident.1 = "L", ident.2 = "H", slot = "counts", test.use = "DESeq2")
de_markers$gene <- rownames(de_markers)
ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,"")), colour = "red", size = 3)
write.table(as.matrix(de_markers),"results/DEG//DESeq2_Seurat_FindMarkers_MUC6LvsH.txt", sep ="\t", col.names = T,row.names = F)

# DEseq2 from matrix, make rank for fgsea----
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
ggplot(res, aes(log2FoldChange, -log10(pvalue))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(padj < 0.1, gene_name,"")), colour = "red", size = 3)
write.table(as.matrix(res),"results/DEG//DESeq2_MUC6LvsH.txt", sep ="\t", col.names = T,row.names = F)

res2 <- res %>% 
  dplyr::select(gene_name, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene_name) %>% 
  summarize(stat=mean(stat))
ranks <- deframe(res2)

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
write.table(fgseaRes[,-8],"results/DEG/fgseaRes_C6_MUC6LvsH.txt", sep ="\t", col.names = T,row.names = F)

plotEnrichment(collections$C6[["EGFR_UP.V1_UP"]], ranks) + labs(title="C6_EGFR_UP.V1_UP")
plotEnrichment(collections$C6[["MEK_UP.V1_UP"]], ranks) + labs(title="C6_MEK_UP.V1_UP")
plotEnrichment(collections$C6[["PDGF_ERK_DN.V1_DN"]], ranks) + labs(title="C6_PDGF_ERK_DN.V1_DN")
plotEnrichment(collections$HALLMARKS[["HALLMARK_KRAS_SIGNALING_DN"]], ranks) + labs(title="HALLMARK_KRAS_SIGNALING_DN")
plotEnrichment(collections$HALLMARKS[["HALLMARK_KRAS_SIGNALING_UP"]], ranks) + labs(title="HALLMARK_KRAS_SIGNALING_UP")

MUC6KO_UP_DN <- gmtPathways("gene_set/MUC6KO_UP_DN_hs.gmt")
fgseaRes <- fgsea(pathways = MUC6KO_UP_DN, stats = ranks, eps=0.0, minSize = 10, maxSize = 500)
plotEnrichment(MUC6KO_UP_DN[["MUC6KO_UP"]], ranks) + labs(title="MUC6KO_UP")
plotEnrichment(MUC6KO_UP_DN[["MUC6KO_DN"]], ranks) + labs(title="MUC6KO_DN")

write.table(fgseaRes[,-8],"results/DEG/fgseaRes_MUC6KO_UP_DN_MUC6LvsH.txt", sep ="\t", col.names = T,row.names = F)
