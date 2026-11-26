rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
while (grDevices::dev.cur() > 1) grDevices::dev.off()
options(stringsAsFactors = FALSE)

SEURAT_RDS <- project_path("results_root", "single_cell/annotation/seurat_object_celltype_0.6.rds")
stopifnot(file.exists(SEURAT_RDS))

workdir <- project_path("results_root", "single_cell/slingshot")
dir.create(workdir, showWarnings = FALSE, recursive = TRUE)
setwd(workdir)

out_prefix <- "preHTC_HTC"

PLOT_W   <- 7
PLOT_H   <- 6
PLOT_DPI <- 300

out_dir_plots  <- normalizePath(file.path(workdir, "plots"),  winslash = "\\", mustWork = FALSE)
out_dir_tables <- normalizePath(file.path(workdir, "tables"), winslash = "\\", mustWork = FALSE)
out_dir_rds    <- normalizePath(file.path(workdir, "rds"),    winslash = "\\", mustWork = FALSE)

dir.create(out_dir_plots,  showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_tables, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_rds,    showWarnings = FALSE, recursive = TRUE)

.plot_path  <- function(fname) file.path(out_dir_plots,  fname)
.table_path <- function(fname) file.path(out_dir_tables, fname)
.rds_path   <- function(fname) file.path(out_dir_rds,    fname)

assay_use <- "RNA"
CELLTYPE_COL <- "cell_type"
SUBSET_CELLTYPES <- c("preHTC", "HTC")

cluster_col <- "cell_type"

SEURAT_UMAP_REDUCTION_FOR_SLING <- "umap"

root_cluster <- "preHTC"
end_clusters <- "HTC"

FILTER_BY_UMAP <- FALSE
FILTER_MODE <- "manual"

UMAP1_MIN <- NA
UMAP1_MAX <- NA
UMAP2_MIN <- NA
UMAP2_MAX <- NA

Q_LOW  <- 0.01
Q_HIGH <- 0.99

PROTECT_ROOT_END_CELLS <- TRUE

DROP_SMALL_CLUSTERS <- TRUE
MIN_CELLS_PER_CLUSTER <- 30
MIN_TOTAL_CELLS_FOR_SLINGSHOT <- 100

ADD_DIRECTION_ARROWS <- TRUE
ARROW_WEIGHT_THRESH  <- 0.60
ARROW_NBINS          <- 14
ARROW_COL            <- "black"
ARROW_LWD            <- 2
ARROW_HEAD_LENGTH    <- 0.08

k_range <- 3:8
evaluate_nGenes <- 500
USE_HVG_ONLY <- TRUE
HVG_N <- 2000
NCORES <- 8

set.seed(1234)

cran_pkgs <- c("ggplot2", "data.table", "Matrix", "matrixStats", "RColorBrewer", "igraph", "FNN")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_pkgs <- c("SingleCellExperiment", "slingshot", "tradeSeq", "SummarizedExperiment", "S4Vectors", "BiocParallel")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(tradeSeq)
  library(SummarizedExperiment)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(matrixStats)
  library(RColorBrewer)
  library(grid)
  library(BiocParallel)
})

.save_plot <- function(p, file_png, w = PLOT_W, h = PLOT_H, dpi = PLOT_DPI) {
  ggsave(filename = file_png, plot = p, width = w, height = h, dpi = dpi, units = "in")
}

.make_pt_cols_blue_red <- function(pt_vec, n = 200) {
  pal <- grDevices::colorRampPalette(c("#0b49a2", "#49a20b", "#a20b49"))(n)
  col <- rep("grey85", length(pt_vec))
  ok <- which(is.finite(pt_vec))
  if (length(ok) < 10) return(col)

  v <- pt_vec[ok]
  br <- unique(as.numeric(stats::quantile(v, probs = seq(0, 1, length.out = n + 1), na.rm = TRUE)))
  if (length(br) < 3) {
    rng <- range(v, na.rm = TRUE)
    br <- seq(rng[1], rng[2], length.out = n + 1)
  }
  idx <- as.integer(cut(v, breaks = br, include.lowest = TRUE, labels = FALSE))
  idx <- pmax(1, pmin(n, idx))
  col[ok] <- pal[idx]
  col
}

.build_arrow_segments_binned <- function(umap_xy, pt_vec, w_vec = NULL,
                                         weight_thresh = 0.6, nbins = 12) {
  stopifnot(nrow(umap_xy) == length(pt_vec))
  ok <- which(is.finite(pt_vec))
  if (!is.null(w_vec)) ok <- intersect(ok, which(is.finite(w_vec) & w_vec >= weight_thresh))
  if (length(ok) < 20) return(list(seg = NULL, lab = NULL, cent = NULL))

  pt_ok <- as.numeric(pt_vec[ok])
  um_ok <- as.matrix(umap_xy[ok, , drop = FALSE])

  qs <- stats::quantile(pt_ok, probs = seq(0, 1, length.out = nbins + 1), na.rm = TRUE)
  breaks <- unique(as.numeric(qs))
  if (length(breaks) < 3) return(list(seg = NULL, lab = NULL, cent = NULL))

  bin_id <- cut(pt_ok, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  bins <- sort(unique(bin_id))
  if (length(bins) < 2) return(list(seg = NULL, lab = NULL, cent = NULL))

  cent <- do.call(rbind, lapply(bins, function(b) {
    idx <- which(bin_id == b)
    data.frame(
      x = mean(um_ok[idx, 1]),
      y = mean(um_ok[idx, 2]),
      bin = b,
      pt_med = stats::median(pt_ok[idx]),
      n = length(idx),
      stringsAsFactors = FALSE
    )
  }))
  cent <- cent[order(cent$pt_med), , drop = FALSE]
  if (nrow(cent) < 2) return(list(seg = NULL, lab = NULL, cent = NULL))

  seg <- data.frame(
    x = cent$x[-nrow(cent)], y = cent$y[-nrow(cent)],
    xend = cent$x[-1], yend = cent$y[-1],
    stringsAsFactors = FALSE
  )
  lab <- data.frame(
    x = c(cent$x[1], cent$x[nrow(cent)]),
    y = c(cent$y[1], cent$y[nrow(cent)]),
    lab = c("early", "late"),
    stringsAsFactors = FALSE
  )
  list(seg = seg, lab = lab, cent = cent)
}

.draw_segments_base <- function(seg_df, lab_df = NULL,
                                arrow_col = "black", arrow_lwd = 2, arrow_head_len = 0.08) {
  if (is.null(seg_df) || nrow(seg_df) == 0) return(invisible(NULL))
  for (i in seq_len(nrow(seg_df))) {
    arrows(seg_df$x[i], seg_df$y[i], seg_df$xend[i], seg_df$yend[i],
           length = arrow_head_len, angle = 25, col = arrow_col, lwd = arrow_lwd)
  }
  if (!is.null(lab_df) && nrow(lab_df) >= 2) {
    text(lab_df$x[1], lab_df$y[1], labels = lab_df$lab[1], pos = 4, cex = 0.9, col = arrow_col)
    text(lab_df$x[2], lab_df$y[2], labels = lab_df$lab[2], pos = 4, cex = 0.9, col = arrow_col)
  }
  invisible(NULL)
}

.get_counts_from_seurat <- function(obj, assay = DefaultAssay(obj)) {
  mat <- tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) NULL
  )
  if (is.null(mat)) {
    mat <- tryCatch(
      Seurat::GetAssayData(obj, assay = assay, slot = "counts"),
      error = function(e) NULL
    )
  }
  if (is.null(mat)) stop("无法从 Seurat 对象中提取 counts 矩阵，请检查 assay/layer。")
  mat
}

.collapse_or_auto <- function(x) {
  if (is.null(x)) return("Auto")
  paste(x, collapse = ",")
}

cat("\n------------------------------\n[STEP] Read Seurat & subset preHTC / HTC\n------------------------------\n")

seurat_full <- readRDS(SEURAT_RDS)
stopifnot(inherits(seurat_full, "Seurat"))

if (!assay_use %in% names(seurat_full@assays)) {
  stop("Seurat 对象缺少 assay=", assay_use, "\n现有 assays: ", paste(names(seurat_full@assays), collapse = ", "))
}
DefaultAssay(seurat_full) <- assay_use

if (!CELLTYPE_COL %in% colnames(seurat_full@meta.data)) {
  stop("Seurat meta.data 中不存在列: ", CELLTYPE_COL)
}
if (!SEURAT_UMAP_REDUCTION_FOR_SLING %in% names(seurat_full@reductions)) {
  stop("Seurat 对象中不存在 reduction: ", SEURAT_UMAP_REDUCTION_FOR_SLING)
}

seurat_full@meta.data[[CELLTYPE_COL]] <- as.character(seurat_full@meta.data[[CELLTYPE_COL]])

cells_subset <- rownames(seurat_full@meta.data)[seurat_full@meta.data[[CELLTYPE_COL]] %in% SUBSET_CELLTYPES]
if (length(cells_subset) == 0) {
  stop("未找到 cell_type 为 ", paste(SUBSET_CELLTYPES, collapse = ", "), " 的细胞。")
}

seurat_sub <- subset(seurat_full, cells = cells_subset)
seurat_sub@meta.data[[CELLTYPE_COL]] <- as.character(seurat_sub@meta.data[[CELLTYPE_COL]])

if (!cluster_col %in% colnames(seurat_sub@meta.data)) {
  stop("subset 后 Seurat meta.data 中不存在 cluster_col: ", cluster_col)
}
seurat_sub@meta.data[[cluster_col]] <- as.character(seurat_sub@meta.data[[cluster_col]])

saveRDS(seurat_sub, .rds_path(paste0(out_prefix, "_subset_preHTC_HTC.rds")))

tab_subset <- data.table(cell_type = seurat_sub@meta.data[[CELLTYPE_COL]])[, .N, by = cell_type][order(-N)]
fwrite(tab_subset, .table_path(paste0(out_prefix, "_00_subset_celltype_counts.csv")))

cat("\n------------------------------\n[STEP] Preview DimPlot & cluster UMAP stats\n------------------------------\n")

p_dim1 <- DimPlot(
  seurat_sub,
  reduction = SEURAT_UMAP_REDUCTION_FOR_SLING,
  group.by = CELLTYPE_COL,
  label = TRUE
) + ggtitle("Subset UMAP by cell_type")
.save_plot(p_dim1, .plot_path(paste0(out_prefix, "_00_DimPlot_umap_celltype.png")), w = 8, h = 6)

p_dim2 <- DimPlot(
  seurat_sub,
  reduction = SEURAT_UMAP_REDUCTION_FOR_SLING,
  group.by = cluster_col,
  label = TRUE
) + ggtitle(paste0("Subset UMAP by ", cluster_col))
.save_plot(p_dim2, .plot_path(paste0(out_prefix, "_00_DimPlot_umap_clustercol.png")), w = 8, h = 6)

umap_full <- Embeddings(seurat_sub, SEURAT_UMAP_REDUCTION_FOR_SLING)[, 1:2, drop = FALSE]
colnames(umap_full) <- c("UMAP_1", "UMAP_2")
cells_all <- rownames(umap_full)

dt_umap <- data.table(
  cell = cells_all,
  cluster = seurat_sub@meta.data[cells_all, cluster_col, drop = TRUE],
  cell_type = seurat_sub@meta.data[cells_all, CELLTYPE_COL, drop = TRUE],
  UMAP_1 = umap_full[cells_all, "UMAP_1"],
  UMAP_2 = umap_full[cells_all, "UMAP_2"]
)
dt_umap <- dt_umap[is.finite(UMAP_1) & is.finite(UMAP_2) & !is.na(cluster) & nzchar(cluster)]

tab_range <- dt_umap[, .(
  n = .N,
  UMAP1_min = min(UMAP_1), UMAP1_max = max(UMAP_1),
  UMAP2_min = min(UMAP_2), UMAP2_max = max(UMAP_2),
  UMAP2_q05 = as.numeric(quantile(UMAP_2, 0.05)),
  UMAP2_q50 = as.numeric(quantile(UMAP_2, 0.50)),
  UMAP2_q95 = as.numeric(quantile(UMAP_2, 0.95))
), by = cluster][order(-n)]
fwrite(tab_range, .table_path(paste0(out_prefix, "_00_cluster_umap_range.csv")))

cat("\n------------------------------\n[STEP] Build fit-cells\n------------------------------\n")

keep_bounds <- rep(TRUE, nrow(dt_umap))
UMAP1_MIN_USE <- NA; UMAP1_MAX_USE <- NA; UMAP2_MIN_USE <- NA; UMAP2_MAX_USE <- NA

if (FILTER_BY_UMAP) {
  if (FILTER_MODE == "quantile") {
    bx <- quantile(dt_umap$UMAP_1, probs = c(Q_LOW, Q_HIGH), na.rm = TRUE)
    by <- quantile(dt_umap$UMAP_2, probs = c(Q_LOW, Q_HIGH), na.rm = TRUE)
    UMAP1_MIN_USE <- as.numeric(bx[1]); UMAP1_MAX_USE <- as.numeric(bx[2])
    UMAP2_MIN_USE <- as.numeric(by[1]); UMAP2_MAX_USE <- as.numeric(by[2])
  } else {
    UMAP1_MIN_USE <- UMAP1_MIN; UMAP1_MAX_USE <- UMAP1_MAX
    UMAP2_MIN_USE <- UMAP2_MIN; UMAP2_MAX_USE <- UMAP2_MAX
  }

  if (!is.na(UMAP1_MIN_USE)) keep_bounds <- keep_bounds & (dt_umap$UMAP_1 >= UMAP1_MIN_USE)
  if (!is.na(UMAP1_MAX_USE)) keep_bounds <- keep_bounds & (dt_umap$UMAP_1 <= UMAP1_MAX_USE)
  if (!is.na(UMAP2_MIN_USE)) keep_bounds <- keep_bounds & (dt_umap$UMAP_2 >= UMAP2_MIN_USE)
  if (!is.na(UMAP2_MAX_USE)) keep_bounds <- keep_bounds & (dt_umap$UMAP_2 <= UMAP2_MAX_USE)

  cat(sprintf("[INFO] Bounds used: UMAP_1 [%s,%s], UMAP_2 [%s,%s]\n",
              as.character(UMAP1_MIN_USE), as.character(UMAP1_MAX_USE),
              as.character(UMAP2_MIN_USE), as.character(UMAP2_MAX_USE)))
}

cells_keep <- dt_umap$cell[keep_bounds]

if (PROTECT_ROOT_END_CELLS) {
  cells_root <- character(0)
  cells_end  <- character(0)
  if (!is.null(root_cluster)) cells_root <- dt_umap$cell[dt_umap$cluster == root_cluster]
  if (!is.null(end_clusters)) cells_end  <- dt_umap$cell[dt_umap$cluster %in% end_clusters]
  cells_keep <- union(cells_keep, union(cells_root, cells_end))
  cat(sprintf("[INFO] PROTECT_ROOT_END_CELLS=TRUE => force keep root/end cells (root=%d, end=%d)\n",
              length(cells_root), length(cells_end)))
}
cells_keep <- unique(cells_keep)

png(.plot_path(paste0(out_prefix, "_01_filter_preview.png")),
    width = PLOT_W, height = PLOT_H, units = "in", res = PLOT_DPI)
plot(dt_umap$UMAP_1, dt_umap$UMAP_2,
     pch = 16, cex = 0.5,
     col = ifelse(dt_umap$cell %in% cells_keep, "black", "grey85"),
     xlab = "UMAP_1", ylab = "UMAP_2",
     main = "Filter preview (black=fit-cells)")
if (FILTER_BY_UMAP && FILTER_MODE != "quantile") {
  if (!is.na(UMAP1_MIN)) abline(v = UMAP1_MIN, lty = 2)
  if (!is.na(UMAP1_MAX)) abline(v = UMAP1_MAX, lty = 2)
  if (!is.na(UMAP2_MIN)) abline(h = UMAP2_MIN, lty = 2)
  if (!is.na(UMAP2_MAX)) abline(h = UMAP2_MAX, lty = 2)
}
dev.off()

tab_after <- data.table(cluster = dt_umap$cluster[dt_umap$cell %in% cells_keep])[, .N, by = cluster][order(-N)]
fwrite(tab_after, .table_path(paste0(out_prefix, "_01_cluster_counts_after_filter.csv")))

if (!is.null(root_cluster)) {
  n_root_after <- tab_after[cluster == root_cluster, N]
  if (length(n_root_after) == 0) n_root_after <- 0
  if (n_root_after < MIN_CELLS_PER_CLUSTER) {
    stop("Root cluster 过滤后细胞过少: ", root_cluster)
  }
}
if (!is.null(end_clusters)) {
  n_end_after <- tab_after[cluster %in% end_clusters, .(cluster, N)]
  if (nrow(n_end_after) == 0 || any(n_end_after$N < MIN_CELLS_PER_CLUSTER)) {
    stop("End cluster 过滤后细胞过少: ", paste(end_clusters, collapse = ","))
  }
}

cat("\n------------------------------\n[STEP] Subset Seurat to fit-cells\n------------------------------\n")

seurat_fit <- subset(seurat_sub, cells = cells_keep)
seurat_fit@meta.data[[cluster_col]] <- as.character(seurat_fit@meta.data[[cluster_col]])

tab_fit <- data.table(cluster = seurat_fit@meta.data[[cluster_col]])[, .N, by = cluster][order(-N)]
fwrite(tab_fit, .table_path(paste0(out_prefix, "_02_cluster_counts_fitcells_beforeDrop.csv")))

if (DROP_SMALL_CLUSTERS) {
  small <- tab_fit[N < MIN_CELLS_PER_CLUSTER, cluster]
  small <- setdiff(small, c(root_cluster, end_clusters))
  if (length(small) > 0) {
    keep2 <- rownames(seurat_fit@meta.data)[!seurat_fit@meta.data[[cluster_col]] %in% small]
    seurat_fit <- subset(seurat_fit, cells = keep2)
    seurat_fit@meta.data[[cluster_col]] <- as.character(seurat_fit@meta.data[[cluster_col]])
  }
}

tab_fit2 <- data.table(cluster = seurat_fit@meta.data[[cluster_col]])[, .N, by = cluster][order(-N)]
fwrite(tab_fit2, .table_path(paste0(out_prefix, "_02_cluster_counts_fitcells_afterDrop.csv")))

if (ncol(seurat_fit) < MIN_TOTAL_CELLS_FOR_SLINGSHOT) {
  stop("Total cells too few for Slingshot: ", ncol(seurat_fit))
}

saveRDS(seurat_fit, .rds_path(paste0(out_prefix, "_Seurat_fitcells_beforeSlingshot.rds")))

cat("\n------------------------------\n[STEP] Ensure HVGs exist\n------------------------------\n")

if (USE_HVG_ONLY) {
  cur_hvg <- tryCatch(VariableFeatures(seurat_fit), error = function(e) character(0))
  if (length(cur_hvg) < 200) {
    cat("[INFO] VariableFeatures too few, running FindVariableFeatures() on subset object.\n")
    seurat_fit <- FindVariableFeatures(
      seurat_fit,
      assay = assay_use,
      selection.method = "vst",
      nfeatures = max(HVG_N, 2000)
    )
  }
}

cat("\n------------------------------\n[STEP] Seurat -> SCE\n------------------------------\n")

counts_seu <- .get_counts_from_seurat(seurat_fit, assay = assay_use)

sce <- as.SingleCellExperiment(seurat_fit, assay = assay_use)

if (!"counts" %in% assayNames(sce)) {
  assay(sce, "counts") <- counts_seu[rownames(sce), colnames(sce), drop = FALSE]
}

umap_fit <- Embeddings(seurat_fit, SEURAT_UMAP_REDUCTION_FOR_SLING)[, 1:2, drop = FALSE]
colnames(umap_fit) <- c("UMAP_1", "UMAP_2")
reducedDim(sce, "UMAP") <- umap_fit
colData(sce)[[cluster_col]] <- droplevels(as.factor(seurat_fit@meta.data[[cluster_col]]))

cat("\n------------------------------\n[STEP] Run Slingshot\n------------------------------\n")

sce_sling <- slingshot(
  sce,
  clusterLabels = cluster_col,
  reducedDim    = "UMAP",
  start.clus    = root_cluster,
  end.clus      = end_clusters
)

SDS <- SlingshotDataSet(sce_sling)
lineage_list <- slingLineages(SDS)
pst_mat <- slingPseudotime(SDS)
wgt_mat <- slingCurveWeights(SDS)

if (ncol(pst_mat) == 0) stop("Slingshot lineage = 0，未得到轨迹。")

topology_file <- .table_path(paste0(out_prefix, "_03_lineage_topology.txt"))
con <- file(topology_file, "wt")
writeLines(sprintf("fit_reduction: UMAP (Seurat %s)", SEURAT_UMAP_REDUCTION_FOR_SLING), con)
writeLines(sprintf("subset cell types: %s", paste(SUBSET_CELLTYPES, collapse = ", ")), con)
writeLines(sprintf("cluster_col  : %s", cluster_col), con)
writeLines(sprintf("root_cluster : %s", .collapse_or_auto(root_cluster)), con)
writeLines(sprintf("end_clusters : %s", .collapse_or_auto(end_clusters)), con)
writeLines("", con)
writeLines("Lineages:", con)
for (i in seq_along(lineage_list)) {
  writeLines(sprintf("Lineage_%d: %s", i, paste(lineage_list[[i]], collapse = " -> ")), con)
}
close(con)

pst_df <- as.data.frame(pst_mat); pst_df$cell <- rownames(pst_df)
wgt_df <- as.data.frame(wgt_mat); wgt_df$cell <- rownames(wgt_df)
fwrite(as.data.table(pst_df), .table_path(paste0(out_prefix, "_03_pseudotime.csv")))
fwrite(as.data.table(wgt_df), .table_path(paste0(out_prefix, "_03_weights.csv")))

n_lin <- ncol(pst_mat)
for (i in seq_len(n_lin)) {
  seurat_fit@meta.data[[paste0("pseudotime_L", i)]] <- pst_mat[, i]
  seurat_fit@meta.data[[paste0("weight_L", i)]]     <- wgt_mat[, i]
}

saveRDS(seurat_fit, .rds_path(paste0(out_prefix, "_Seurat_with_slingshot.rds")))
saveRDS(sce_sling,  .rds_path(paste0(out_prefix, "_SCE_with_slingshot.rds")))

cat("\n------------------------------\n[STEP] Plot (Arrows Only, No Curves) for ALL lineages\n------------------------------\n")

umap_xy <- reducedDim(sce_sling, "UMAP")[, 1:2, drop = FALSE]
colnames(umap_xy) <- c("UMAP_1", "UMAP_2")

fname_base <- paste0(out_prefix, "_UMAPmain_FULL_ArrowsOnly_AllLineages.png")
png(.plot_path(fname_base), width = PLOT_W, height = PLOT_H, units = "in", res = PLOT_DPI)
plot(umap_xy[, 1], umap_xy[, 2], pch = 16, cex = 0.55, asp = 1,
     xlab = "UMAP_1", ylab = "UMAP_2", col = "grey80",
     main = "Slingshot (Arrows Only) - All Lineages")

if (ADD_DIRECTION_ARROWS) {
  for (i in seq_len(n_lin)) {
    pt <- pst_mat[, i]
    w  <- wgt_mat[, i]
    arrow_obj <- .build_arrow_segments_binned(
      umap_xy = umap_xy,
      pt_vec = pt,
      w_vec  = w,
      weight_thresh = ARROW_WEIGHT_THRESH,
      nbins = ARROW_NBINS
    )
    if (!is.null(arrow_obj$seg)) {
      .draw_segments_base(arrow_obj$seg, arrow_obj$lab, ARROW_COL, ARROW_LWD, ARROW_HEAD_LENGTH)
    }
  }
}
dev.off()

for (i in seq_len(n_lin)) {

  pt_i <- pst_mat[, i]
  w_i  <- wgt_mat[, i]
  pt_col <- .make_pt_cols_blue_red(pt_i, n = 200)

  arrow_obj <- .build_arrow_segments_binned(
    umap_xy = umap_xy,
    pt_vec = pt_i,
    w_vec  = w_i,
    weight_thresh = ARROW_WEIGHT_THRESH,
    nbins = ARROW_NBINS
  )

  fname_li <- paste0(out_prefix, "_UMAPmain_FULL_pseudotime_L", i, ".png")

  png(.plot_path(fname_li), width = PLOT_W, height = PLOT_H, units = "in", res = PLOT_DPI)
  plot(umap_xy[, 1], umap_xy[, 2], pch = 16, cex = 0.55, asp = 1,
       col = pt_col,
       xlab = "UMAP_1", ylab = "UMAP_2",
       main = paste0("Lineage ", i, ": Pseudotime (Arrows Only)"))

  if (ADD_DIRECTION_ARROWS && !is.null(arrow_obj$seg)) {
    .draw_segments_base(arrow_obj$seg, arrow_obj$lab, ARROW_COL, ARROW_LWD, ARROW_HEAD_LENGTH)
  }

  dev.off()
}

cat("\n------------------------------\n[STEP] tradeSeq\n------------------------------\n")

if (!"counts" %in% assayNames(sce_sling)) {
  stop("SCE 中没有 assay 'counts'。")
}
counts <- assay(sce_sling, "counts")

if (USE_HVG_ONLY) {
  hvgs <- NULL
  cur_hvg <- tryCatch(VariableFeatures(seurat_fit), error = function(e) character(0))

  if (length(cur_hvg) > 0) {
    hvgs <- intersect(cur_hvg, rownames(counts))
    hvgs <- hvgs[seq_len(min(HVG_N, length(hvgs)))]
  }

  if (is.null(hvgs) || length(hvgs) < 200) {
    rv <- matrixStats::rowVars(as.matrix(counts))
    keep_genes <- order(rv, decreasing = TRUE)[seq_len(min(HVG_N, length(rv)))]
    counts_use <- counts[keep_genes, , drop = FALSE]
  } else {
    counts_use <- counts[hvgs, , drop = FALSE]
  }
} else {
  counts_use <- counts
}

icMat <- evaluateK(
  counts  = counts_use,
  sds     = SDS,
  k       = k_range,
  nGenes  = min(evaluate_nGenes, nrow(counts_use)),
  verbose = TRUE
)

ic_mean <- colMeans(icMat, na.rm = TRUE)
k_from_col <- suppressWarnings(as.integer(gsub("[^0-9]", "", colnames(icMat))))
if (length(k_from_col) != length(ic_mean) || all(is.na(k_from_col))) {
  k_from_col <- as.integer(k_range)
}
k_from_col <- k_from_col[seq_len(min(length(k_from_col), length(ic_mean)))]
ic_mean <- ic_mean[seq_len(length(k_from_col))]

best_k <- k_from_col[which.min(ic_mean)]
if (length(best_k) == 0 || is.na(best_k)) {
  best_k <- as.integer(stats::median(k_range))
}

fwrite(
  data.table(k = k_from_col, IC_mean = ic_mean, best = (k_from_col == best_k)),
  .table_path(paste0(out_prefix, "_tradeSeq_evaluateK_summary.csv"))
)

pseudotime  <- slingPseudotime(SDS, na = FALSE)
cellWeights <- slingCurveWeights(SDS)

fitGAM_formals <- names(formals(tradeSeq::fitGAM))
fitGAM_args <- list(
  counts      = counts_use,
  pseudotime  = pseudotime,
  cellWeights = cellWeights,
  nknots      = best_k,
  verbose     = TRUE
)

if ("ncores" %in% fitGAM_formals) {
  fitGAM_args$ncores <- NCORES
}
if ("BPPARAM" %in% fitGAM_formals) {
  if (.Platform$OS.type == "windows") {
    fitGAM_args$BPPARAM <- BiocParallel::SnowParam(workers = max(1, NCORES), type = "SOCK")
  } else {
    fitGAM_args$BPPARAM <- BiocParallel::MulticoreParam(workers = max(1, NCORES))
  }
}

sce_gam <- do.call(tradeSeq::fitGAM, fitGAM_args)
saveRDS(sce_gam, .rds_path(paste0(out_prefix, "_tradeSeq_fitGAM.rds")))

assoc_res <- associationTest(sce_gam)
assoc_res$gene <- rownames(assoc_res)

sve_res <- startVsEndTest(sce_gam)
sve_res$gene <- rownames(sve_res)

fwrite(as.data.table(assoc_res), .table_path(paste0(out_prefix, "_tradeSeq_associationTest.csv")))
fwrite(as.data.table(sve_res),   .table_path(paste0(out_prefix, "_tradeSeq_startVsEndTest.csv")))

if (n_lin >= 2) {
  tr_ver <- tryCatch(as.character(packageVersion("tradeSeq")), error = function(e) NA_character_)
  fwrite(data.table(tradeSeq_version = tr_ver),
         .table_path(paste0(out_prefix, "_tradeSeq_version.csv")))

  has_branchTest  <- FALSE
  has_patternTest <- FALSE
  has_diffEndTest <- FALSE

  if (requireNamespace("tradeSeq", quietly = TRUE)) {
    exps <- tryCatch(getNamespaceExports("tradeSeq"), error = function(e) character(0))
    has_branchTest  <- "branchTest"  %in% exps
    has_patternTest <- "patternTest" %in% exps
    has_diffEndTest <- "diffEndTest" %in% exps
  }

  if (has_branchTest) {
    message("[INFO] tradeSeq::branchTest is available -> run branchTest()")
    br_res <- tradeSeq::branchTest(sce_gam)
    br_res$gene <- rownames(br_res)
    fwrite(as.data.table(br_res), .table_path(paste0(out_prefix, "_tradeSeq_branchTest.csv")))
  } else {
    message("[WARN] tradeSeq::branchTest NOT found. Use alternatives.")
    if (has_patternTest) {
      pt_res <- tradeSeq::patternTest(sce_gam)
      pt_res$gene <- rownames(pt_res)
      fwrite(as.data.table(pt_res), .table_path(paste0(out_prefix, "_tradeSeq_patternTest.csv")))
    }
    if (has_diffEndTest) {
      de_res <- tradeSeq::diffEndTest(sce_gam)
      de_res$gene <- rownames(de_res)
      fwrite(as.data.table(de_res), .table_path(paste0(out_prefix, "_tradeSeq_diffEndTest.csv")))
    }
  }
}

cat("\n------------------------------\n[STEP] Write unified pseudotime back to Seurat\n------------------------------\n")

LINEAGE_TO_PLOT <- 1
seurat_fit@meta.data[["pseudotime"]] <- seurat_fit@meta.data[[paste0("pseudotime_L", LINEAGE_TO_PLOT)]]
saveRDS(seurat_fit, .rds_path(paste0(out_prefix, "_Seurat_with_pseudotimeUnified.rds")))

p_fp <- FeaturePlot(
  seurat_fit,
  features = "pseudotime",
  reduction = SEURAT_UMAP_REDUCTION_FOR_SLING
)
.save_plot(p_fp, .plot_path(paste0(out_prefix, "_Seurat_FeaturePlot_pseudotime_on_umap.png")))

cat("============================================================\n")
cat("Slingshot + tradeSeq pipeline DONE (Arrows Only Mode).\n")
cat(sprintf("Input Seurat : %s\n", SEURAT_RDS))
cat(sprintf("Subset types : %s\n", paste(SUBSET_CELLTYPES, collapse = ", ")))
cat(sprintf("cluster_col  : %s\n", cluster_col))
cat(sprintf("root_cluster : %s\n", .collapse_or_auto(root_cluster)))
cat(sprintf("end_clusters : %s\n", .collapse_or_auto(end_clusters)))
cat(sprintf("Lineages     : %d\n", n_lin))
cat(sprintf("Workdir      : %s\n", workdir))
cat("============================================================\n")
