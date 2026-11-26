rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
while (grDevices::dev.cur() > 1) grDevices::dev.off()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

SEURAT_RDS <- project_path("data_root", "auxiliary/trajectory/annotated_seurat.rds")
WORKDIR    <- project_path("results_root", "cytotrace")
stopifnot(file.exists(SEURAT_RDS))
dir.create(WORKDIR, showWarnings = FALSE, recursive = TRUE)
setwd(WORKDIR)

SPECIES   <- "human"
ASSAY_USE <- "RNA"
SLOT_USE  <- "counts"
SEED      <- 14

NCORES            <- 1
BATCH_SIZE        <- 5000
SMOOTH_BATCH_SIZE <- 1000

GROUP_COL <- "cell_type"
REDUCTION_PRIORITY <- c("umap_main", "umap")

seu <- readRDS(SEURAT_RDS)
DefaultAssay(seu) <- ASSAY_USE

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

if (!requireNamespace("CytoTRACE2", quietly = TRUE)) {
  remotes::install_github(
    repo = "digitalcytometry/cytotrace2",
    subdir = "cytotrace2_r",
    dependencies = TRUE,
    upgrade = "never",
    build_vignettes = FALSE
  )
}
suppressPackageStartupMessages(library(CytoTRACE2))

if (!ASSAY_USE %in% names(seu@assays)) {
  stop(paste0("[ERROR] assay 不存在：", ASSAY_USE, "；可选：", paste(names(seu@assays), collapse = ", ")))
}
if (!SLOT_USE %in% c("counts", "data", "scale.data")) {
  stop("[ERROR] SLOT_USE 仅允许 counts / data / scale.data")
}

expr0 <- GetAssayData(seu, assay = ASSAY_USE, slot = SLOT_USE)

if (inherits(expr0, "dgCMatrix")) {

  if (length(expr0@x) > 0 && any(expr0@x < 0)) stop("[ERROR] 检测到负值：CytoTRACE2 不允许负值。请用 raw counts/CPM/TPM（非log）")

  q99 <- stats::quantile(expr0@x, probs = 0.99, na.rm = TRUE)
  if (!is.na(q99) && q99 < 20) {
    message("[WARN] 你的矩阵 99% 分位 < 20，看起来可能是 log/data slot。CytoTRACE2 推荐 counts/raw（非log）。如确实是 raw counts 可忽略。")
  }
} else {

  if (any(expr0 < 0, na.rm = TRUE)) stop("[ERROR] 检测到负值：CytoTRACE2 不允许负值。请用 raw counts/CPM/TPM（非log）")
  q99 <- stats::quantile(as.numeric(expr0), probs = 0.99, na.rm = TRUE)
  if (!is.na(q99) && q99 < 20) {
    message("[WARN] 你的矩阵 99% 分位 < 20，看起来可能是 log/data slot。CytoTRACE2 推荐 counts/raw（非log）。如确实是 raw counts 可忽略。")
  }
}

seu_ct2 <- CytoTRACE2::cytotrace2(
  input                 = seu,
  species               = SPECIES,
  is_seurat             = TRUE,
  slot_type             = SLOT_USE,
  batch_size            = BATCH_SIZE,
  smooth_batch_size     = SMOOTH_BATCH_SIZE,
  parallelize_models    = (NCORES > 1),
  parallelize_smoothing = (NCORES > 1),
  ncores                = NCORES,
  seed                  = SEED
)

reducs <- names(seu_ct2@reductions)
if (length(reducs) == 0) stop("[ERROR] Seurat 对象里没有任何 reductions，无法画 UMAP/Embedding 图。")

reduction_use <- NULL
for (r in REDUCTION_PRIORITY) {
  if (r %in% reducs) { reduction_use <- r; break }
}
if (is.null(reduction_use)) reduction_use <- reducs[1]
message("[INFO] 使用降维：", reduction_use)

if (!"CytoTRACE2_Score" %in% colnames(seu_ct2@meta.data)) {
  stop("[ERROR] meta.data 中未发现 CytoTRACE2_Score；cytotrace2() 可能未成功写入。")
}

p_score <- FeaturePlot(
  object = seu_ct2,
  features = "CytoTRACE2_Score",
  reduction = reduction_use
) + ggtitle("CytoTRACE2_Score")

ggsave(
  filename = file.path(WORKDIR, "UMAP_CytoTRACE2_Score.png"),
  plot = p_score, width = 7, height = 6, dpi = 300
)

if ("CytoTRACE2_Potency" %in% colnames(seu_ct2@meta.data)) {
  p_pot <- DimPlot(
    object = seu_ct2,
    group.by = "CytoTRACE2_Potency",
    reduction = reduction_use
  ) + ggtitle("CytoTRACE2_Potency")

  ggsave(
    filename = file.path(WORKDIR, "UMAP_CytoTRACE2_Potency.png"),
    plot = p_pot, width = 7, height = 6, dpi = 300
  )
} else {
  message("[WARN] meta.data 中未发现 CytoTRACE2_Potency（部分版本/参数可能不输出该列），跳过 potency 图。")
}

if (GROUP_COL %in% colnames(seu_ct2@meta.data)) {
  p_vln <- VlnPlot(
    object = seu_ct2,
    features = "CytoTRACE2_Score",
    group.by = GROUP_COL,
    pt.size = 0
  ) + ggtitle(paste0("CytoTRACE2_Score by ", GROUP_COL))

  ggsave(
    filename = file.path(WORKDIR, paste0("Vln_CytoTRACE2_Score_by_", GROUP_COL, ".png")),
    plot = p_vln, width = 10, height = 5, dpi = 300
  )
} else {
  message("[INFO] meta.data 中没有列：", GROUP_COL, "；跳过小提琴图。")
}

saveRDS(seu_ct2, file = file.path(WORKDIR, "seu_with_CytoTRACE2.rds"))

write.csv(
  x = seu_ct2@meta.data,
  file = file.path(WORKDIR, "meta_with_CytoTRACE2.csv"),
  row.names = TRUE
)

cat("\n[DONE]\n")
cat("Workdir: ", WORKDIR, "\n", sep = "")
cat("Saved:   seu_with_CytoTRACE2.rds\n")
cat("Plots:   UMAP_CytoTRACE2_Score.png, UMAP_CytoTRACE2_Potency.png (if exists), Vln_* (if exists)\n")
