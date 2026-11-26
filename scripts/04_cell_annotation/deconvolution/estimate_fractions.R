rm(list = ls())

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
options(stringsAsFactors = FALSE)

library(Seurat)
library(BisqueRNA)
library(Biobase)
library(dplyr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

seurat_file <- project_path("data_root", "auxiliary/oa/fibrocartilage_annotated.rds")

bulk_count_file <- project_path("data_root", "bulk_deconvolution/OA_COUNTS_raw.csv")

bulk_meta_file  <- project_path("data_root", "bulk_deconvolution/sample_info.csv")

state_col   <- "cell_type2"
patient_col <- "patient"

states_keep <- c("EC", "HomC", "HTC",
                 "MMP13-", "MMP13+",
                 "preHTC", "proC",
                 "RegC", "RepC")

cat(">>> 读取 bulk 计数矩阵...\n")

bulk_raw <- read.csv(
  bulk_count_file,
  header      = TRUE,
  row.names   = 1,
  check.names = FALSE
)

bulk.counts <- as.matrix(bulk_raw)

cat("bulk 计数矩阵维度: ",
    nrow(bulk.counts), " genes × ", ncol(bulk.counts), " samples\n")

cat(">>> 读取 Seurat 对象...\n")
obj <- readRDS(seurat_file)
DefaultAssay(obj) <- "RNA"

meta <- obj@meta.data

if (!state_col %in% colnames(meta)) {
  stop("在 Seurat meta 中找不到列: ", state_col)
}
if (!patient_col %in% colnames(meta)) {
  stop("在 Seurat meta 中找不到列: ", patient_col)
}

cells_use <- rownames(meta)[meta[[state_col]] %in% states_keep]
cat("保留下来的细胞数: ", length(cells_use), "\n")

obj_sub <- subset(obj, cells = cells_use)
rm(obj); gc()

sc_counts_sparse <- GetAssayData(obj_sub, assay = "RNA", slot = "counts")

genes_bulk   <- rownames(bulk.counts)
genes_sc_all <- rownames(sc_counts_sparse)
genes_use    <- intersect(genes_bulk, genes_sc_all)

genes_use <- genes_use[
  Matrix::rowSums(sc_counts_sparse[genes_use, , drop = FALSE] > 0) >= 10
]

cat("最终保留的基因数：", length(genes_use), "\n")

sc.counts <- as.matrix(sc_counts_sparse[genes_use, , drop = FALSE])
cat("单细胞计数矩阵维度: ",
    nrow(sc.counts), " genes × ", ncol(sc.counts), " cells\n")

rm(sc_counts_sparse); gc()

meta_sub <- obj_sub@meta.data
rm(obj_sub)
gc()
stopifnot(all(colnames(sc.counts) == rownames(meta_sub)))

sc.pheno <- data.frame(
  SubjectName = as.character(meta_sub[[patient_col]]),
  cellType    = as.character(meta_sub[[state_col]]),
  row.names   = colnames(sc.counts),
  stringsAsFactors = FALSE
)

sc.varMeta <- data.frame(
  labelDescription = c("individual / patient id", "cell state / type"),
  row.names        = c("SubjectName", "cellType"),
  stringsAsFactors = FALSE
)

sc.pdata <- new("AnnotatedDataFrame",
                data        = sc.pheno,
                varMetadata = sc.varMeta)

sc.eset <- ExpressionSet(
  assayData = sc.counts,
  phenoData = sc.pdata
)
rm(sc.counts)
gc()
cat("sc.eset 构建完成。\n")

bulk.eset <- ExpressionSet(assayData = bulk.counts)

cat(">>> 从 bulk.counts 列名推断样本信息（group / region）...\n")

samples <- colnames(bulk.counts)

bulk.meta <- data.frame(
  sample = samples,
  stringsAsFactors = FALSE
)

bulk.meta$group <- NA_character_
bulk.meta$group[grepl("^Control", bulk.meta$sample)] <- "HC"
bulk.meta$group[grepl("^(NWB|WB)", bulk.meta$sample)] <- "OA"

library(stringr)
region_mat <- str_match(bulk.meta$sample, "\\.(AS|SZ|MZ|DZ)")
bulk.meta$region <- region_mat[, 2]

bulk.meta$group  <- factor(bulk.meta$group,  levels = c("HC", "OA"))
bulk.meta$region <- factor(bulk.meta$region, levels = c("AS", "SZ", "MZ", "DZ"))

cat("bulk 样本信息构建完成：\n")
print(head(bulk.meta))

idx <- match(colnames(bulk.counts), bulk.meta$sample)

bulk.meta <- bulk.meta[idx, , drop = FALSE]

bulk.meta$group  <- factor(bulk.meta$group,  levels = c("HC", "OA"))
bulk.meta$region <- factor(bulk.meta$region, levels = c("AS", "SZ", "MZ", "DZ"))

cat("bulk 样本信息对齐完成。\n")

cat(">>> 使用 BisqueRNA::ReferenceBasedDecomposition 反卷积...\n")

res.bisque <- ReferenceBasedDecomposition(
  bulk.eset,
  sc.eset,
  markers     = NULL,
  use.overlap = FALSE
)

bulk.props <- res.bisque$bulk.props

cat("反卷积完成，得到比例矩阵: ",
    nrow(bulk.props), " samples × ", ncol(bulk.props), " cell states\n")

out_prop_file <- project_path("data_root", "bulk_deconvolution/Bisque_cell_fraction_matrix.csv")
write.csv(bulk.props, file = out_prop_file)
cat("细胞比例矩阵已保存到: ", out_prop_file, "\n")

rm(list = ls())
options(stringsAsFactors = FALSE)

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

prop_file <- project_path("data_root", "bulk_deconvolution/Bisque_cell_fraction_matrix.csv")

cat(">>> 读取 Bisque 结果矩阵...\n")
prop_mat <- read.csv(
  prop_file,
  row.names   = 1,
  check.names = FALSE
)

cat("原始矩阵维度: ", nrow(prop_mat), " × ", ncol(prop_mat), "\n")

bulk.props <- t(as.matrix(prop_mat))

cat("转置后维度（样本 × 细胞类型）: ",
    nrow(bulk.props), " samples × ", ncol(bulk.props), " cell types\n")

row_sum_head <- head(rowSums(bulk.props))
print(row_sum_head)

props_df <- as.data.frame(bulk.props)
props_df$sample <- rownames(props_df)

props_df$group <- ifelse(grepl("^Control", props_df$sample),
                         "HC", "OA")
props_df$group <- factor(props_df$group, levels = c("HC", "OA"))

table(props_df$group)

props_long <- props_df |>
  pivot_longer(
    cols      = -c(sample, group),
    names_to  = "cell_type",
    values_to = "prop"
  )

summary(props_long$prop)

props_summary <- props_long |>
  group_by(group, cell_type) |>
  summarise(
    prop = mean(prop, na.rm = TRUE),
    .groups = "drop"
  )

celltype_levels <- c("EC", "HomC", "HTC",
                     "MMP13-", "MMP13+",
                     "preHTC", "proC",
                     "RegC", "RepC")

props_summary$cell_type <- factor(
  props_summary$cell_type,
  levels = intersect(celltype_levels, unique(props_summary$cell_type))
)

cat(">>> 全流程结束。\n")
