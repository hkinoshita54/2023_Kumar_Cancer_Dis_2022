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

# PB group by clusters (cell_type)----
seu_epi_t <- readRDS("RDSfiles/seu_epi_t.RDS")
Idents(seu_epi_t) <- "cell_type"
DimPlot(seu_epi_t)
bulk <- bulk <- AggregateExpression(seu_epi_t, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("cell_type", "sample_ID", "MUC6_level"))
tail(Cells(bulk))
bulk$cell_type <- sapply(strsplit(Cells(bulk), split = "_"), "[", 1)
bulk$sample_ID <- sapply(strsplit(Cells(bulk), split = "_"), "[", 2)
bulk$MUC6_level <- sapply(strsplit(Cells(bulk), split = "-"), "[", 2)
bulk$MUC6_level <- factor(x = bulk$MUC6_level, levels = c("L", "H"))

# Intestinal----
Intestinal.bulk <- subset(bulk, cell_type == "Intestinal")
Idents(Intestinal.bulk) <- "MUC6_level"
de_markers <- FindMarkers(Intestinal.bulk, ident.1 = "L", ident.2 = "H", slot = "counts", test.use = "DESeq2", verbose = F)
de_markers$gene <- rownames(de_markers)
ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.1, gene, "")), colour = "red", size = 3)

# DESeq2 from count matrix----
cts <- as.matrix(Intestinal.bulk[["RNA"]]$counts)
conditions <- sapply(strsplit(Cells(Intestinal.bulk), split = "-"), "[", 2)
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
# write.table(as.matrix(res),"./Results/DEG//DESeq2_MUC6LvsH_1007.txt", sep ="\t", col.names = T,row.names = F)
ggplot(res, aes(log2FoldChange, -log10(pvalue))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(padj < 0.01, gene_name,"")), colour = "red", size = 3)

res2 <- res %>% 
  dplyr::select(gene_name, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(gene_name) %>% 
  summarize(stat=mean(stat))
ranks <- deframe(res2)

# fgsea----
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
