suppressPackageStartupMessages({
  library(Seurat)
  library(Nebulosa)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

if (!all(c("group", "weight") %in% colnames(seurat_object@meta.data))) {
  stop("meta.data 里必须包含列: 'group' 和 'weight'。")
}

cells_ALL    <- WhichCells(seurat_object)
cells_HC     <- WhichCells(seurat_object, expression = group == "HC")
cells_OA     <- WhichCells(seurat_object, expression = group == "OA")
cells_OA_WB  <- WhichCells(seurat_object, expression = group == "OA" & weight == "WB")
cells_OA_NWB <- WhichCells(seurat_object, expression = group == "OA" & weight == "NWB")

plot_density_fixed <- function(obj, cells, title, feature = "MMP13",
                               reduction_use = "umap",
                               max_den = 0.35) {
  if (length(cells) == 0) {
    message("状态 [", title, "] 没有细胞，跳过。")
    return(NULL)
  }

  obj_sub <- subset(obj, cells = cells)

  p <- Nebulosa::plot_density(
    obj_sub,
    features  = feature,
    reduction = reduction_use
  ) +
    scale_fill_viridis_c(
      limits = c(0, max_den),
      oob    = scales::squish
    ) +
    scale_colour_viridis_c(
      limits = c(0, max_den),
      oob    = scales::squish
    ) +
    ggtitle(title) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )

  return(p)
}

plots_list <- list(
  ALL    = plot_density_fixed(seurat_object, cells_ALL,     "ALL - MMP13"),
  HC     = plot_density_fixed(seurat_object, cells_HC,     "HC - MMP13"),
  OA     = plot_density_fixed(seurat_object, cells_OA,     "OA - MMP13"),
  OA_WB  = plot_density_fixed(seurat_object, cells_OA_WB,  "OA WB - MMP13"),
  OA_NWB = plot_density_fixed(seurat_object, cells_OA_NWB, "OA NWB - MMP13")
)

plots_list <- Filter(Negate(is.null), plots_list)

if (length(plots_list) == 0) {
  stop("四个状态都没有可用细胞，无法作图。")
}

p_all_mmp13 <- wrap_plots(plots_list, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

print(p_all_mmp13)

suppressPackageStartupMessages({
  library(Seurat)
  library(Nebulosa)
  library(ggplot2)
  library(patchwork)
})

if (!all(c("group", "weight") %in% colnames(seurat_object@meta.data))) {
  stop("meta.data 里必须包含列: 'group' 和 'weight'。")
}

cells_ALL    <- WhichCells(seurat_object)
cells_HC     <- WhichCells(seurat_object, expression = group == "HC")
cells_OA     <- WhichCells(seurat_object, expression = group == "OA")
cells_OA_WB  <- WhichCells(seurat_object, expression = group == "OA" & weight == "WB")
cells_OA_NWB <- WhichCells(seurat_object, expression = group == "OA" & weight == "NWB")

plot_density_auto <- function(obj, cells, title,
                              feature = "MMP13",
                              reduction_use = "umap") {
  if (length(cells) == 0) {
    message("状态 [", title, "] 没有细胞，跳过。")
    return(NULL)
  }

  obj_sub <- subset(obj, cells = cells)

  p <- Nebulosa::plot_density(
    obj_sub,
    features  = feature,
    reduction = reduction_use
  ) +
    ggtitle(title) +
    theme(
      plot.title      = element_text(hjust = 0.5),
      legend.position = "right"
    )

  return(p)
}

plots_list <- list(
  ALL    = plot_density_auto(seurat_object, cells_ALL,     "ALL - MMP13"),
  HC     = plot_density_auto(seurat_object, cells_HC,      "HC - MMP13"),
  OA     = plot_density_auto(seurat_object, cells_OA,      "OA - MMP13"),
  OA_WB  = plot_density_auto(seurat_object, cells_OA_WB,   "OA WB - MMP13"),
  OA_NWB = plot_density_auto(seurat_object, cells_OA_NWB,  "OA NWB - MMP13")
)

plots_list <- Filter(Negate(is.null), plots_list)

if (length(plots_list) == 0) {
  stop("四个状态都没有可用细胞，无法作图。")
}

p_all_mmp13 <- wrap_plots(plots_list, ncol = 2)

print(p_all_mmp13)
