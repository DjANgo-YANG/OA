rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
while (grDevices::dev.cur() > 1) grDevices::dev.off()
options(stringsAsFactors = FALSE)

SEURAT_WITH_SLING_RDS <- project_path("data_root", "auxiliary/trajectory/slingshot/rds/FINAL_LIG_Slingshot_tradeSeq_umapMainFit_root11_end6_FilteredRobust_Seurat_with_slingshot.rds")

workdir <- project_path("data_root", "auxiliary/trajectory/slingshot")
out_dir_tables <- normalizePath(file.path(workdir, "tables"), winslash = "\\", mustWork = FALSE)
out_dir_rds    <- normalizePath(file.path(workdir, "rds"),    winslash = "\\", mustWork = FALSE)
dir.create(out_dir_tables, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_rds,    showWarnings = FALSE, recursive = TRUE)

prefix <- "FINAL_LIG_cluster6_signature_for_NicheNet"

.table_path <- function(fname) file.path(out_dir_tables, fname)
.rds_path   <- function(fname) file.path(out_dir_rds, fname)

assay_use   <- "RNA"
cluster_col <- "seurat_clusters"

receiver_cluster <- "6"
root_cluster     <- "11"

min_pct         <- 0.10
logfc_threshold <- 0.25
only_pos        <- TRUE
test_use        <- "wilcox"

top_n_geneset <- 50

pct_in_cluster6 <- 0.10

model_dir <- project_path("references", "nichenet_v2")
organism  <- "human"

set.seed(1)

need_cran <- c("data.table", "Matrix", "dplyr", "tibble")
for (p in need_cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tibble)
})

stopifnot(file.exists(SEURAT_WITH_SLING_RDS))
seu <- readRDS(SEURAT_WITH_SLING_RDS)
stopifnot(inherits(seu, "Seurat"))

if (!assay_use %in% names(seu@assays)) stop("缺少 assay=", assay_use)
DefaultAssay(seu) <- assay_use

if (!cluster_col %in% colnames(seu@meta.data)) {
  stop("meta.data 缺少 cluster_col=", cluster_col,
       "\n可用列: ", paste(colnames(seu@meta.data), collapse = ", "))
}
seu@meta.data[[cluster_col]] <- as.character(seu@meta.data[[cluster_col]])

cts <- sort(unique(seu@meta.data[[cluster_col]]))
if (!receiver_cluster %in% cts) stop("receiver_cluster 不存在：", receiver_cluster)
if (!root_cluster %in% cts) stop("root_cluster 不存在：", root_cluster)

.read_nichenet_ltm <- function(model_dir, organism = c("human", "mouse")) {
  organism <- match.arg(organism)
  if (!dir.exists(model_dir)) stop("model_dir 不存在：", model_dir)

  if (organism == "human") {
    p_ltm <- file.path(model_dir, "ligand_target_matrix_nsga2r_final.rds")
  } else {
    p_ltm <- file.path(model_dir, "ligand_target_matrix_nsga2r_final_mouse.rds")
  }
  if (!file.exists(p_ltm)) stop("找不到 ligand_target_matrix：", p_ltm)

  ltm <- readRDS(p_ltm)
  if (!is.matrix(ltm) && !inherits(ltm, "Matrix")) stop("ligand_target_matrix 结构不对")
  if (is.null(rownames(ltm))) stop("ligand_target_matrix 缺 rownames (targets)")
  ltm
}

ligand_target_matrix <- .read_nichenet_ltm(model_dir, organism)
targets_space <- rownames(ligand_target_matrix)

.get_expressed_pct <- function(seu, cells, assay = NULL, slot = "data", pct = 0.10) {
  if (is.null(assay)) assay <- DefaultAssay(seu)
  expr <- GetAssayData(seu, assay = assay, slot = slot)[, cells, drop = FALSE]
  det <- Matrix::rowSums(expr != 0)
  pct_vec <- det / length(cells)
  names(pct_vec)[pct_vec >= pct]
}

cells_c6 <- rownames(seu@meta.data)[seu@meta.data[[cluster_col]] == receiver_cluster]
expr_genes_c6 <- .get_expressed_pct(seu, cells_c6, assay = assay_use, slot = "data", pct = pct_in_cluster6)

Idents(seu) <- seu@meta.data[[cluster_col]]

mk_c6_vs_all <- FindMarkers(
  object = seu,
  ident.1 = receiver_cluster,
  ident.2 = NULL,
  assay = assay_use,
  slot = "data",
  test.use = test_use,
  min.pct = min_pct,
  logfc.threshold = logfc_threshold,
  only.pos = only_pos
) %>% rownames_to_column("gene")

fc_col <- c("avg_log2FC", "avg_logFC")[c("avg_log2FC", "avg_logFC") %in% colnames(mk_c6_vs_all)][1]
if (is.na(fc_col)) stop("FindMarkers 结果没有 avg_log2FC/avg_logFC 列")

mk_c6_vs_all <- mk_c6_vs_all %>%
  mutate(logFC = suppressWarnings(as.numeric(.data[[fc_col]]))) %>%
  arrange(p_val_adj, desc(logFC))

geneset_c6_vs_all_raw <- mk_c6_vs_all$gene

mk_c6_vs_11 <- FindMarkers(
  object = seu,
  ident.1 = receiver_cluster,
  ident.2 = root_cluster,
  assay = assay_use,
  slot = "data",
  test.use = test_use,
  min.pct = min_pct,
  logfc.threshold = logfc_threshold,
  only.pos = only_pos
) %>% rownames_to_column("gene")

fc_col2 <- c("avg_log2FC", "avg_logFC")[c("avg_log2FC", "avg_logFC") %in% colnames(mk_c6_vs_11)][1]
if (is.na(fc_col2)) stop("FindMarkers(c6 vs 11) 结果没有 avg_log2FC/avg_logFC 列")

mk_c6_vs_11 <- mk_c6_vs_11 %>%
  mutate(logFC = suppressWarnings(as.numeric(.data[[fc_col2]]))) %>%
  arrange(p_val_adj, desc(logFC))

geneset_c6_vs_11_raw <- mk_c6_vs_11$gene

.filter_geneset <- function(genes, expr_genes, targets_space, top_n = 300) {
  genes <- unique(genes)
  genes <- intersect(genes, expr_genes)
  genes <- intersect(genes, targets_space)
  head(genes, top_n)
}

geneset_c6_vs_all <- .filter_geneset(geneset_c6_vs_all_raw, expr_genes_c6, targets_space, top_n_geneset)
geneset_c6_vs_11  <- .filter_geneset(geneset_c6_vs_11_raw,  expr_genes_c6, targets_space, top_n_geneset)

fwrite(as.data.table(mk_c6_vs_all), .table_path(paste0(prefix, "_marker_cluster6_vs_all.csv")))
fwrite(as.data.table(mk_c6_vs_11),  .table_path(paste0(prefix, "_marker_cluster6_vs_11.csv")))

fwrite(data.table(gene = geneset_c6_vs_all),
       .table_path(paste0(prefix, "_geneset_cluster6_vs_all.csv")))
fwrite(data.table(gene = geneset_c6_vs_11),
       .table_path(paste0(prefix, "_geneset_cluster6_vs_11.csv")))

meta_out <- data.table(
  seurat_rds = SEURAT_WITH_SLING_RDS,
  assay_use = assay_use,
  cluster_col = cluster_col,
  receiver_cluster = receiver_cluster,
  root_cluster = root_cluster,
  min_pct = min_pct,
  logfc_threshold = logfc_threshold,
  only_pos = only_pos,
  pct_in_cluster6 = pct_in_cluster6,
  top_n_geneset = top_n_geneset,
  geneset_cluster6_vs_all_n = length(geneset_c6_vs_all),
  geneset_cluster6_vs_11_n  = length(geneset_c6_vs_11)
)
fwrite(meta_out, .table_path(paste0(prefix, "_geneset_cluster6_meta.csv")))

saveRDS(
  list(
    params = as.list(meta_out[1]),
    geneset_cluster6_vs_all = geneset_c6_vs_all,
    geneset_cluster6_vs_11  = geneset_c6_vs_11,
    expressed_genes_cluster6 = expr_genes_c6,
    marker_table_cluster6_vs_all = mk_c6_vs_all,
    marker_table_cluster6_vs_11  = mk_c6_vs_11
  ),
  .rds_path(paste0(prefix, "_geneset_for_nichenet.rds"))
)

cat("============================================================\n")
cat("[DONE] NicheNet geneset files written:\n")
cat("  ", .table_path(paste0(prefix, "_geneset_cluster6_vs_all.csv")), "\n", sep="")
cat("  ", .table_path(paste0(prefix, "_geneset_cluster6_vs_11.csv")), "\n", sep="")
cat("Choose ONE of them as geneset_oi in NicheNet.\n")
cat("============================================================\n")
