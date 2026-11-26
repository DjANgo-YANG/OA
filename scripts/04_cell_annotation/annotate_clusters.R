rm(list=ls())
gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))

library(SingleR)
library(celldex)
library(Seurat)
library(BiocGenerics)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(scrapper)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(wesanderson)
library(dplyr)
library(Nebulosa)
library(SCP)

work_dir <- project_path("results_root", "single_cell/annotation")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)
seurat_object <- readRDS(project_path("results_root", "single_cell/qc/4.rds"))
Idents(seurat_object) <- "RNA_snn_res.0.6"
SEURAT_OBJECT <- seurat_object$seurat_clusters <- seurat_object$RNA_snn_res.0.6

all.markers <- FindAllMarkers(
  seurat_object,
  only.pos = TRUE,
  min.pct = 0.4,
  logfc.threshold = 0.25
)
write.csv(all.markers,file = "ALLMARKERS0.6.csv")

head(all.markers)

top5 <- all.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)
write.csv(top5,file = "top5_0.6.csv")

top.genes <- unique(top5$gene)
top.genes
seurat_object@meta.data$seurat_clusters <- seurat_object$RNA_snn_res.0.6

Idents(seurat_object) <- "seurat_clusters"
p <- DoHeatmap(
  seurat_object,
  features = top.genes,
  group.by = "seurat_clusters",
  size =3,
  label = TRUE
) +
  theme(axis.text.y = element_text(size = 6))

p

pdf("heatmap_top5_0.6.pdf",width = 7,height = 5)
p
dev.off

library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)
library(rlang)

marker_groups <- list(
  RegC_1 = c("MT1G", "MT1H"),
  RegC_2 = c("CHI3L2", "MT1X"),
  EC = c("CYTL1", 'FRZB', 'CHAD', 'CLEC3A', 'S100A1', 'S100B'),
  preFC = c('PRG4', 'SEMA3A', 'OGN', 'COL3A1'),
  FC1 = c('MMP2','IGFBP5', 'S100A4', 'MSMP', 'COL1A1', 'COL1A2'),
  FC2 = c('MMP2','IGFBP5', 'S100A4', 'PRG4', 'SEMA3A'),
  preHTC = c('SOX9', 'COL9A3',"COL2A1"),
  HTC = c('COL10A1', 'IBSP', 'SOX9', 'JUN','SPP1'),
  HomC = c('MMP3', 'RGS16','HSPA1B', 'HSPA1A'),
  RepC = c('CILP', 'CILP2', 'COL2A1')
)

panel_order <- c(
  "RegC_1","EC","preHTC","FC2","preFC","RegC_2","HomC","HTC","FC1","RepC"
)

marker_groups <- marker_groups[panel_order]
feat_ordered <- unlist(marker_groups, use.names = FALSE)
grp <- rep(names(marker_groups), lengths(marker_groups))
gene2grp <- data.frame(
  features.plot = feat_ordered,
  panel = factor(grp, levels = panel_order),
  stringsAsFactors = FALSE
)

unique_feats <- unique(feat_ordered)

p <- DotPlot(
  seurat_object,
  features  = unique_feats,
  dot.scale = 8,
  group.by  = "cell_type"
)

p$data <- left_join(p$data, gene2grp, by = "features.plot")
p$data <- dplyr::filter(p$data, !is.na(panel))

p$data <- p$data %>%
  mutate(feature_panel = paste0(panel, "___", features.plot))

feature_panel_levels <- with(gene2grp, paste0(panel, "___", features.plot))

feature_panel_levels <- unique(feature_panel_levels)
p$data$feature_panel <- factor(p$data$feature_panel, levels = feature_panel_levels)

labels <- setNames(gsub("^.*___", "", feature_panel_levels), feature_panel_levels)

p2 <- ggplot(p$data, aes(x = feature_panel, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  facet_grid(~ panel, scales = "free_x", space = "free_x") +
  scale_color_gradient2(
    low       = "#0808ff",
    mid       = "#ffffff",
    high      = "#fd2220",
    midpoint  = 0,
    limits    = c(-2, 2),
    oob       = squish,
    breaks    = c(-2, 0, 2),
    name      = "Average Expression"
  ) +

  scale_size(range = c(0, 6), name = "Pct Expressed") +
  scale_x_discrete(labels = labels) +
  theme_bw() +
  theme(
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.spacing.x    = unit(1.2, "lines"),
    strip.background   = element_rect(fill = "grey92", color = "black"),
    strip.text.x       = element_text(face = "bold",
                                      margin = margin(t = 2, r = 4, b = 2, l = 4)),
    axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = element_line(color = "grey95")
  ) +
  labs(x = NULL, y = NULL)
p2
pdf("DotPlot_marker_genes.pdf",width = 20,height = 6)
p2
dev.off()

seurat_object@meta.data$cell_type <- recode(seurat_object$seurat_clusters,
                                            "0" = "RegC_1",
                                            "1" = "EC",
                                            "2" = "preHTC",
                                            "3" = "FC2",
                                            "4" = "preFC",
                                            "5" = "RegC_2",
                                            "6" = "HomC",
                                            "7" = "HTC",
                                            "8" = "FC1",
                                            "9" = "RepC",
                                            "10" = "RegC_2"

)
table(seurat_object@meta.data$cell_type)

p2 = CellDimPlot(seurat_object, group.by = "cell_type", reduction = "UMAP",
                 pt.size = 0.8)
pdf("UMAP_celltype_0.6.pdf",width = 7,height = 5)
p2
dev.off()
plot_density(seurat_object,features = "MMP13",reduction = "umap",size = 0.3)

p3 = VlnPlot(seurat_object, features = "ADAMTS5", group.by = "cell_type",pt.size = 0)
p3
saveRDS(seurat_object,file = "seurat_object_celltype_0.6.rds")
