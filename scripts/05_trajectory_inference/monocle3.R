rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
while (grDevices::dev.cur() > 1) grDevices::dev.off()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(SeuratWrappers)
  library(ggplot2)
  library(data.table)
})

SEURAT_RDS <- project_path("data_root", "auxiliary/oa/mmp13_subset.rds")
stopifnot(file.exists(SEURAT_RDS))

OUTDIR <- project_path("results_root", "monocle3")
PREFIX <- "MMP13_traj"

CLUSTER_COL <- "cell_type"

ROOT_GROUP_NAME <- "MMP13-"
ROOT_MAX_CELLS  <- 2000

PRUNE_GRAPH <- FALSE
MIN_BRANCH_LEN <- 5

PLOT_W <- 8
PLOT_H <- 6
PLOT_DPI <- 300

set.seed(1234)

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
out_plots  <- file.path(OUTDIR, "plots");  dir.create(out_plots,  showWarnings = FALSE, recursive = TRUE)
out_rds    <- file.path(OUTDIR, "rds");    dir.create(out_rds,    showWarnings = FALSE, recursive = TRUE)
out_tables <- file.path(OUTDIR, "tables"); dir.create(out_tables, showWarnings = FALSE, recursive = TRUE)

.plot_path  <- function(fname) file.path(out_plots,  fname)
.rds_path   <- function(fname) file.path(out_rds,    fname)
.table_path <- function(fname) file.path(out_tables, fname)

message("Reading Seurat object: ", SEURAT_RDS)
seu <- readRDS(SEURAT_RDS)

if (!CLUSTER_COL %in% colnames(seu@meta.data)) {
  stop("Seurat meta.data 中找不到列：", CLUSTER_COL, "\n现有列：", paste(colnames(seu@meta.data), collapse=", "))
}

message("Converting Seurat to Monocle3 cell_data_set...")
cds <- as.cell_data_set(seu)

if ("umap" %in% names(seu@reductions)) {
  reducedDims(cds)$UMAP <- Embeddings(seu, "umap")
}

recreate_partition <- function(g) {

  rep("1", length(g))
}

tryCatch({
  cds <- cluster_cells(cds, reduction_method = "UMAP")
}, error = function(e) {
  message("cluster_cells 失败，尝试手动初始化 partitions...")
  partitions(cds) <- as.factor(rep("1", ncol(cds)))
  names(partitions(cds)) <- colnames(cds)
})

meta <- colData(cds)
root_candidates <- rownames(meta)[meta[[CLUSTER_COL]] == ROOT_GROUP_NAME]

if (length(root_candidates) == 0) {
  stop("在 ", CLUSTER_COL, " 列中找不到名为 '", ROOT_GROUP_NAME, "' 的细胞！\n",
       "该列的取值有：", paste(unique(meta[[CLUSTER_COL]]), collapse=", "))
}

message("Found ", length(root_candidates), " root cells (", ROOT_GROUP_NAME, ")")

root_cells <- root_candidates
if (length(root_cells) > ROOT_MAX_CELLS) {
  root_cells <- sample(root_cells, ROOT_MAX_CELLS)
}

p_diag <- plot_cells(cds, color_cells_by = CLUSTER_COL, show_trajectory_graph = FALSE) +
  ggtitle(paste0("Diagnostic: ", CLUSTER_COL))
ggsave(.plot_path("00_Diagnostic_Cluster.png"), p_diag, width = PLOT_W, height = PLOT_H, dpi = PLOT_DPI)

message("Running FIX_A: learn_graph(use_partition = FALSE)...")

cds_A <- cds
cds_A <- learn_graph(
  cds_A,
  use_partition = FALSE,
  learn_graph_control = list(
    prune_graph = PRUNE_GRAPH,
    minimal_branch_len = MIN_BRANCH_LEN
  )
)

cds_A <- order_cells(cds_A, root_cells = root_cells)

pA1 <- plot_cells(cds_A, color_cells_by = CLUSTER_COL,
                  label_groups_by_cluster = FALSE,
                  label_leaves = FALSE, label_branch_points = FALSE) +
  ggtitle("Trajectory (Color by CellType)")
ggsave(.plot_path("FIX_A_01_CellType.png"), pA1, width = PLOT_W, height = PLOT_H, dpi = PLOT_DPI)

pA2 <- plot_cells(cds_A, color_cells_by = "pseudotime",
                  label_cell_groups = FALSE,
                  label_leaves = FALSE, label_branch_points = FALSE) +
  ggtitle(paste0("Pseudotime (Root: ", ROOT_GROUP_NAME, ")"))
ggsave(.plot_path("FIX_A_02_Pseudotime.png"), pA2, width = PLOT_W, height = PLOT_H, dpi = PLOT_DPI)

saveRDS(cds_A, .rds_path("FIX_A_cds.rds"))

pt_df <- data.frame(
  cell_id = colnames(cds_A),
  cell_type = colData(cds_A)[[CLUSTER_COL]],
  pseudotime = pseudotime(cds_A)
)
write.csv(pt_df, .table_path("FIX_A_pseudotime_values.csv"), row.names = FALSE)

target_group <- "MMP13+"
if (target_group %in% unique(pt_df$cell_type)) {

  p_box <- ggplot(pt_df, aes(x = cell_type, y = pseudotime, fill = cell_type)) +
    geom_boxplot() +
    theme_classic() +
    ggtitle("Pseudotime Distribution by Group") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(.plot_path("FIX_A_03_Pseudotime_Boxplot.png"), p_box, width = 6, height = 5)

  message("已生成起点与终点的 Pseudotime 箱线图对比。")
}

message("========================================")
message("Done! 结果保存在：", OUTDIR)
message("重点查看 plots/FIX_A_02_Pseudotime.png 确认轨迹方向。")
