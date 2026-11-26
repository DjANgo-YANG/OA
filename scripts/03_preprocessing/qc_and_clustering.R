rm(list=ls())
gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))

library(dplyr)
library(tidyr)
library(patchwork)
library(Seurat)
library(ggplot2)
work_dir <- project_path("results_root", "single_cell/qc")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)
seurat_object <- readRDS(project_path("results_root", "single_cell/doublets/nodouble.rds"))
Idents(seurat_object) <- seurat_object@meta.data$orig.ident

seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
seurat_object[["percent.ribo"]] <- PercentageFeatureSet(seurat_object, pattern = "^RP[SL]")
seurat_object[["percent.hb"]] <- PercentageFeatureSet(seurat_object, pattern = "^HB[^(P)]")

pdf("nFeature_RNA-nCount_RNA.pdf", width = 7, height = 5)
VlnPlot(seurat_object,
        features = c("nFeature_RNA", "nCount_RNA"),
        ncol = 2,
        pt.size = 0)
dev.off()

pdf("nFeature_RNA-nCount_RNA-percent.mt.pdf", width = 7, height = 5)
VlnPlot(seurat_object,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        pt.size = 0)
dev.off()

pdf("percent.ribo-percent.hb.pdf", width = 7, height = 5)
VlnPlot(seurat_object,
        features = c("percent.ribo", "percent.hb"),
        ncol = 2,
        pt.size = 0)
dev.off()

obj <- seurat_object

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

obj$All <- "scRNA"
Idents(obj) <- "scRNA"

pdf("ALL-nFeature_RNA-nCount_RNA-percent.mt.pdf", width = 7, height = 5)
VlnPlot(
  obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "All",
  pt.size = 0,
  cols = "#eb6a6a"
) + NoLegend()
dev.off()

plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.ribo")

pdf("相关性分析.pdf", width = 15, height = 5)
plot1 + plot2+plot3
dev.off()

seurat_object <- subset(seurat_object, subset = nFeature_RNA > 500 & nFeature_RNA < 6400 & nCount_RNA < 40000 &  percent.mt < 8 & percent.ribo < 29 )

quantile(seurat_object$nCount_RNA,seq(0.01,0.1,0.01))
quantile(seurat_object$nCount_RNA,seq(0.9,1,0.01))
quantile(seurat_object$nFeature_RNA,seq(0.01,0.1,0.01))
quantile(seurat_object$nFeature_RNA,seq(0.9,1,0.01))
quantile(seurat_object$percent.mt,seq(0.9,1,0.01))
quantile(seurat_object$percent.ribo,seq(0.9,1,0.01))

pdf("nFeature_RNA-nCount_RNA-percent.mt-percent.ribo-percent.hb.pdf", width = 14, height = 10)
VlnPlot(seurat_object,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo", "percent.hb"),
        pt.size = 0)
dev.off()

obj2 <- seurat_object

obj2[["percent.mt"]] <- PercentageFeatureSet(obj2, pattern = "^MT-")

obj2$All <- "scRNA"
Idents(obj2) <- "scRNA"

pdf("all-nFeature_RNA-nCount_RNA-percent.mt.pdf", width = 14, height = 5)
VlnPlot(
  obj2,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "All",
  pt.size = 0,
  cols = "#eb6a6a"
) + NoLegend()
dev.off()

plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.ribo")

pdf("相关性分析.pdf", width = 14, height = 5)
plot1 + plot2+plot3
dev.off()
saveRDS(seurat_object,file ="3.rds")

seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(seurat_object), 10)
plot1 <- VariableFeaturePlot(seurat_object)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

combined_plot <- plot1 + plot2 + plot_layout(ncol = 2, guides = 'collect') & theme(legend.position = 'bottom')
pdf("前10个高变基因.pdf", width = 14, height = 5)
combined_plot
dev.off()
print(combined_plot)

CaseMatch(c(cc.genes$s.genes,cc.genes$g2m.genes),VariableFeatures(seurat_object))
g2m_genes = CaseMatch(search = cc.genes$g2m.genes,match = rownames(seurat_object))
s_genes = CaseMatch(search = cc.genes$s.genes,match = rownames(seurat_object))

seurat_object <- CellCycleScoring(seurat_object,s.features = s_genes,g2m.features = g2m_genes)
table(seurat_object$Phase)

seurat_object[["RNA"]] <- split(seurat_object[["RNA"]], f = seurat_object$sample)

seurat_object <- FindVariableFeatures(seurat_object, nfeatures=3000)
seurat_object <- ScaleData(seurat_object,
                           vars.to.regress=c("S.Score","G2M.Score","percent.mt"),
                           features=VariableFeatures(seurat_object))
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

pdf("去批次前DIMPLOT.pdf", width = 7, height = 5)
DimPlot(object = seurat_object,reduction = "pca",pt.size = 1,group.by = "sample")
dev.off()

library(harmony)
seurat_object <- RunHarmony(seurat_object,
                            group.by.vars="sample",
                            plot_convergence = TRUE,
                            max.iter = 60,
                            theta = 2)

DimPlot(object = seurat_object,reduction = "harmony",pt.size = 1,group.by = "sample")

pdf("去批次后DIMPLOT.pdf", width = 7, height = 5)
DimPlot(object = seurat_object,reduction = "harmony",pt.size = 1,group.by = "sample")
dev.off()

seurat_object <- JoinLayers(seurat_object)

ElbowPlot(seurat_object,ndims = 50,reduction = "pca")

stdev <- seurat_object[["pca"]]@stdev

var_explained <- stdev^2 / sum(stdev^2)

cum_var <- cumsum(var_explained)

plota <- plot(cum_var,
     xlab = "PC number",
     ylab = "Cumulative Proportion of Variance Explained",
     type = "b", pch = 19, col = "blue")
abline(h = 0.7, col = "red", lty = 2)
abline(h = 0.8, col = "darkgreen", lty = 2)
pdf("累计方差解释率曲线.pdf", width = 7, height = 5)
plota
dev.off()

seurat_object <- FindNeighbors(seurat_object,reduction = "harmony",dims = 1:20)

seurat_object <- FindClusters(seurat_object, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2))
table(seurat_object$RNA_snn_res.0.6)
Idents(seurat_object) <- "RNA_snn_res.0.6"

seurat_object <- RunUMAP(seurat_object,reduction = 'harmony',
                         dims        = 1:20,
                         min.dist = 0.6,
                         n.neighbors = 40,
                         spread = 1)
p1 <- DimPlot(seurat_object, reduction = "umap",label = T)
pdf("UMAP聚类图0.6.pdf", width = 7, height = 5)
p1
dev.off()
saveRDS(seurat_object,file = "4.rds")
print(p1)

p_cc <- DimPlot(
  seurat_object,
  reduction = "umap",
  group.by = "Phase",
  label = FALSE, repel = TRUE,
  raster = TRUE
)

pdf("细胞周期UMAP.pdf", width = 7, height = 5)
p_cc
dev.off()

p_oid <- DimPlot(
  seurat_object,
  reduction = "umap",
  group.by = "orig.ident",
  label = FALSE,
  raster = TRUE
)
pdf("样本来源UMAP.pdf", width = 7, height = 5)
p_oid
dev.off()
