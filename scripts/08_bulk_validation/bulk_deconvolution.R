rm(list = ls())

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(SingleCellExperiment)
  library(S4Vectors)
  library(MuSiC)

  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(stringr)
})

work_dir <- project_path("data_root", "bulk_reference")
setwd(work_dir)

seurat_file     <- project_path("data_root", "auxiliary/reference/annotation/seurat_object_celltype_0.6.rds")
bulk_count_file <- project_path("data_root", "bulk_reference/OA_COUNTS_raw.csv")

out_dir <- file.path(work_dir, "MuSiC_official_mergedTypes_clean_oriented")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

state_col   <- "cell_type"
patient_col <- "orig.ident"
min_cells   <- 30

is_bulk_sample <- function(x) grepl("^(Control|WB|NWB)\\.", x)

cat(">>> Read bulk counts...\n")
bulk_raw <- read.csv(bulk_count_file, header = TRUE, row.names = 1, check.names = FALSE)
bulk.mtx0 <- as.matrix(bulk_raw)

storage.mode(bulk.mtx0) <- "numeric"
bulk.mtx0[is.na(bulk.mtx0) | is.nan(bulk.mtx0) | is.infinite(bulk.mtx0)] <- 0

lib_bulk0 <- colSums(bulk.mtx0)
if (any(lib_bulk0 == 0)) {
  bad <- names(lib_bulk0)[lib_bulk0 == 0]
  cat("Dropping bulk samples with library size 0:\n", paste(bad, collapse = "\n"), "\n")
  bulk.mtx0 <- bulk.mtx0[, lib_bulk0 > 0, drop = FALSE]
}
cat("bulk dim: ", nrow(bulk.mtx0), " genes × ", ncol(bulk.mtx0), " samples\n")

cat(">>> Read Seurat...\n")
obj <- readRDS(seurat_file)
DefaultAssay(obj) <- "RNA"

meta <- obj@meta.data
if (!state_col %in% colnames(meta)) stop("找不到列: ", state_col)
if (!patient_col %in% colnames(meta)) stop("找不到列: ", patient_col)

ct0 <- as.character(meta[[state_col]])
ct1 <- ct0
ct1[ct0 %in% c("RegC_1", "RegC_2")] <- "RegC"
ct1[ct0 %in% c("preHTC", "HTC")]    <- "HT_like"
ct1[ct0 %in% c("preFC", "FC1", "FC2")] <- "FC_like"
meta$cell_type_merged <- ct1
obj@meta.data$cell_type_merged <- ct1

cnt0 <- GetAssayData(obj, assay = "RNA", slot = "counts")

genes_use <- intersect(rownames(cnt0), rownames(bulk.mtx0))
cat("overlap genes:", length(genes_use), "\n")

cnt <- cnt0[genes_use, , drop = FALSE]
bulk.mtx <- bulk.mtx0[genes_use, , drop = FALSE]

sc.sce <- SingleCellExperiment(
  assays  = list(counts = cnt),
  colData = S4Vectors::DataFrame(
    SubjectName = as.character(meta[[patient_col]]),
    cellType    = as.character(meta$cell_type_merged),
    row.names   = colnames(cnt)
  )
)

rm(obj, cnt0, bulk_raw, bulk.mtx0); gc()

ct_tab <- table(colData(sc.sce)$cellType)
cat("cellType counts:\n"); print(ct_tab)

keep_ct <- names(ct_tab)[ct_tab >= min_cells]
if (length(keep_ct) < length(ct_tab)) {
  cat("Dropping rare cell types (<", min_cells, " cells): ",
      paste(setdiff(names(ct_tab), keep_ct), collapse = ", "), "\n")
  sc.sce <- sc.sce[, colData(sc.sce)$cellType %in% keep_ct]
  colData(sc.sce)$cellType <- droplevels(factor(colData(sc.sce)$cellType))
}

sc_counts <- assay(sc.sce, "counts")
keep_gene <- Matrix::rowSums(sc_counts) > 0 & rowSums(bulk.mtx) > 0

bulk.mtx2 <- bulk.mtx[keep_gene, , drop = FALSE]
sc.sce2   <- sc.sce[keep_gene, ]

cat("After cleaning: genes=", nrow(bulk.mtx2),
    " bulk samples=", ncol(bulk.mtx2),
    " sc cells=", ncol(sc.sce2), "\n")

cat(">>> Run MuSiC::music_prop...\n")
res <- MuSiC::music_prop(
  bulk.mtx = bulk.mtx2,
  sc.sce   = sc.sce2,
  clusters = "cellType",
  samples  = "SubjectName",
  verbose  = TRUE
)

prop0 <- res$Est.prop.weighted
cat("MuSiC output dim:", nrow(prop0), "×", ncol(prop0), "\n")
cat("MuSiC rownames head:\n"); print(head(rownames(prop0)))
cat("MuSiC colnames head:\n"); print(head(colnames(prop0)))

rn_hits <- sum(is_bulk_sample(rownames(prop0)))
cn_hits <- sum(is_bulk_sample(colnames(prop0)))
cat(">>> Orientation check: row hits=", rn_hits, " col hits=", cn_hits, "\n")

if (cn_hits > rn_hits) {

  prop_sample_by_ct <- t(prop0)
  prop_ct_by_sample <- prop0
  cat("Using transpose: prop0 is celltype×sample\n")
} else if (rn_hits > cn_hits) {

  prop_sample_by_ct <- prop0
  prop_ct_by_sample <- t(prop0)
  cat("No transpose: prop0 is sample×celltype\n")
} else {

  if (nrow(prop0) == ncol(bulk.mtx2)) {
    prop_sample_by_ct <- prop0
    prop_ct_by_sample <- t(prop0)
    cat("Fallback decided: treat prop0 as sample×celltype\n")
  } else {
    prop_sample_by_ct <- t(prop0)
    prop_ct_by_sample <- prop0
    cat("Fallback decided: treat prop0 as celltype×sample\n")
  }
}

write.csv(prop_ct_by_sample, file.path(out_dir, "MuSiC_prop_celltype_by_sample.csv"))
write.csv(prop_sample_by_ct, file.path(out_dir, "MuSiC_prop_sample_by_celltype.csv"))

plot_df <- as.data.frame(prop_sample_by_ct, check.names = FALSE)
plot_df$sample <- rownames(plot_df)

if (sum(is_bulk_sample(plot_df$sample)) == 0) {
  stop("自动判向后 sample 仍不像 bulk 样本名。请检查 MuSiC 输出的行/列名。")
}

plot_df$group <- ifelse(grepl("^Control\\.", plot_df$sample), "HC", "OA")
plot_df$group <- factor(plot_df$group, levels = c("HC", "OA"))

region_mat <- stringr::str_match(plot_df$sample, "\\.(AS|SZ|MZ|DZ)")
plot_df$region <- region_mat[, 2]
plot_df$region <- factor(plot_df$region, levels = c("AS", "SZ", "MZ", "DZ"))

if (any(is.na(plot_df$region))) {
  bad <- plot_df$sample[is.na(plot_df$region)]
  stop("以下样本无法解析 region(AS/SZ/MZ/DZ)：\n", paste(bad, collapse = "\n"))
}

plot_long <- plot_df %>%
  pivot_longer(cols = -c(sample, group, region),
               names_to = "cell_type",
               values_to = "prop")

cell_levels <- c("EC", "RegC", "FC_like", "HT_like", "HomC", "RepC")
plot_long$cell_type <- factor(plot_long$cell_type, levels = intersect(cell_levels, unique(plot_long$cell_type)))

write.csv(plot_long, file.path(out_dir, "MuSiC_props_per_sample_long.csv"), row.names = FALSE)

p1_df <- plot_long %>%
  group_by(group, cell_type) %>%
  summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop")

p1 <- ggplot(p1_df, aes(x = group, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = percent_format()) +
  theme_bw(base_size = 12) +
  labs(x = "", y = "Cell fraction", fill = "Cell type")

pdf1 <- file.path(out_dir, "HC_vs_OA_stacked.pdf")
pdf(pdf1, width = 6, height = 5)
print(p1)
dev.off()

p2_df <- plot_long %>%
  group_by(group, region, cell_type) %>%
  summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop")

p2 <- ggplot(p2_df, aes(x = group, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ region, nrow = 1) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw(base_size = 12) +
  labs(x = "", y = "Cell fraction", fill = "Cell type")

pdf2 <- file.path(out_dir, "groupXregion_stacked.pdf")
pdf(pdf2, width = 10, height = 4.5)
print(p2)
dev.off()

cat("Saved outputs to:\n", out_dir, "\n>>> DONE.\n")
