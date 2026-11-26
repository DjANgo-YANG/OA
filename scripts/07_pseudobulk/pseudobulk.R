rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
while (grDevices::dev.cur() > 1) grDevices::dev.off()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(DESeq2)
})

SEURAT_RDS <- project_path("data_root", "auxiliary/trajectory/annotated_seurat.rds")
stopifnot(file.exists(SEURAT_RDS))

WORKDIR <- project_path("data_root", "auxiliary/trajectory/pseudobulk")
dir.create(WORKDIR, showWarnings = FALSE, recursive = TRUE)
setwd(WORKDIR)

GENES <- c("TGFB1", "IL1B")

ASSAY_USE <- "RNA"
SLOT_USE  <- "counts"

SAMPLE_COL   <- "auto"
GROUP_COL    <- "auto"
CELLTYPE_COL <- "auto"

DO_CELLTYPE_PSEUDOBULK <- FALSE

GROUP_LEVELS_PREFERRED <- c("PLL", "OPLL")

SIG_CUTS <- c(0.001, 0.01, 0.05)

seu <- readRDS(SEURAT_RDS)
stopifnot(inherits(seu, "Seurat"))

if (!ASSAY_USE %in% names(seu@assays)) {
  message("[WARN] ASSAY_USE 不存在：", ASSAY_USE, "；可用 assays: ", paste(names(seu@assays), collapse = ", "))
  if ("RNA" %in% names(seu@assays)) ASSAY_USE <- "RNA"
  else if ("SCT" %in% names(seu@assays)) ASSAY_USE <- "SCT"
  else ASSAY_USE <- names(seu@assays)[1]
  message("[INFO] 已回退 ASSAY_USE = ", ASSAY_USE)
}

meta <- seu@meta.data

.pick_col <- function(candidates, meta_df) {
  hit <- candidates[candidates %in% colnames(meta_df)]
  if (length(hit) > 0) return(hit[1])
  return(NA_character_)
}

if (identical(SAMPLE_COL, "auto")) {
  SAMPLE_COL <- .pick_col(candidates = c("orig.ident", "sample", "Sample"), meta_df = meta)
}
if (identical(GROUP_COL, "auto")) {
  GROUP_COL <- .pick_col(candidates = c("group", "condition", "status"), meta_df = meta)
}
if (identical(CELLTYPE_COL, "auto")) {
  CELLTYPE_COL <- .pick_col(candidates = c("cell_type", "celltype", "annotation", "seurat_clusters"), meta_df = meta)
}

if (is.na(SAMPLE_COL) || is.na(GROUP_COL)) {
  stop(
    "无法自动识别 SAMPLE_COL 或 GROUP_COL。\n",
    "请在脚本 0.4 手动指定列名：SAMPLE_COL / GROUP_COL。\n",
    "当前 meta 列名如下：\n",
    paste(colnames(meta), collapse = ", ")
  )
}
if (DO_CELLTYPE_PSEUDOBULK && is.na(CELLTYPE_COL)) {
  stop(
    "你开启了 DO_CELLTYPE_PSEUDOBULK=TRUE，但无法识别 CELLTYPE_COL。\n",
    "请在脚本 0.4 手动指定 CELLTYPE_COL。\n",
    "当前 meta 列名如下：\n",
    paste(colnames(meta), collapse = ", ")
  )
}

message("[INFO] SAMPLE_COL = ", SAMPLE_COL)
message("[INFO] GROUP_COL  = ", GROUP_COL)
message("[INFO] CELLTYPE_COL = ", CELLTYPE_COL)

DefaultAssay(seu) <- ASSAY_USE
counts_mat <- GetAssayData(seu, assay = ASSAY_USE, slot = SLOT_USE)

genes_present <- GENES[GENES %in% rownames(counts_mat)]
genes_missing <- setdiff(GENES, genes_present)

if (length(genes_missing) > 0) {
  message("[WARN] 以下基因在表达矩阵中未找到：", paste(genes_missing, collapse = ", "))
  for (g in genes_missing) {
    hits <- grep(g, rownames(counts_mat), ignore.case = TRUE, value = TRUE)
    if (length(hits) > 0) {
      message("  - 可能相关符号 (", g, "): ", paste(head(hits, 20), collapse = ", "),
              ifelse(length(hits) > 20, " ...", ""))
    }
  }
}
if (length(genes_present) == 0) stop("目标基因全部未在矩阵中找到；请检查基因符号（人：TGFB1/IL1B；鼠：Tgfb1/Il1b 等）。")

meta_use <- meta %>%
  mutate(
    .sample   = as.character(.data[[SAMPLE_COL]]),
    .group    = as.character(.data[[GROUP_COL]]),
    .celltype = if (!is.na(CELLTYPE_COL)) as.character(.data[[CELLTYPE_COL]]) else NA_character_
  )

keep_cells <- rownames(meta_use)[!is.na(meta_use$.sample) & meta_use$.sample != "" &
                                   !is.na(meta_use$.group)  & meta_use$.group  != ""]
if (length(keep_cells) < 10) stop("可用细胞数过少（<10）；请检查 meta 里 SAMPLE_COL/GROUP_COL 是否正确。")

seu2 <- subset(seu, cells = keep_cells)
meta2 <- seu2@meta.data %>%
  mutate(
    .sample   = as.character(.data[[SAMPLE_COL]]),
    .group    = as.character(.data[[GROUP_COL]]),
    .celltype = if (!is.na(CELLTYPE_COL)) as.character(.data[[CELLTYPE_COL]]) else NA_character_
  )

if (!DO_CELLTYPE_PSEUDOBULK) {
  meta2$.agg_id <- meta2$.sample
} else {
  meta2$.agg_id <- paste(meta2$.sample, meta2$.celltype, sep = "__")
}
seu2@meta.data$.agg_id <- meta2$.agg_id

pb_list <- AggregateExpression(
  seu2,
  assays = ASSAY_USE,
  slot = SLOT_USE,
  group.by = ".agg_id",
  return.seurat = FALSE
)
pb_counts <- as.matrix(pb_list[[ASSAY_USE]])

agg_map <- meta2 %>%
  dplyr::select(.agg_id, .sample, .group, .celltype) %>%
  distinct()

bad_agg <- agg_map %>%
  group_by(.agg_id) %>%
  summarise(n_group = n_distinct(.group), .groups = "drop") %>%
  filter(n_group > 1)
if (nrow(bad_agg) > 0) {
  stop("同一个聚合单元 .agg_id 对应多个 group（PLL/OPLL 混杂）；请检查 meta 的 GROUP_COL。\n",
       paste(bad_agg$.agg_id, collapse = ", "))
}

coldata <- agg_map %>%
  group_by(.agg_id) %>%
  summarise(
    sample   = .sample[1],
    group    = .group[1],
    celltype = if (DO_CELLTYPE_PSEUDOBULK) .celltype[1] else "ALL",
    .groups  = "drop"
  )

coldata <- coldata %>% filter(.agg_id %in% colnames(pb_counts))
pb_counts <- pb_counts[, coldata$.agg_id, drop = FALSE]

groups_in_data <- unique(coldata$group)
lvl <- GROUP_LEVELS_PREFERRED[GROUP_LEVELS_PREFERRED %in% groups_in_data]
if (length(lvl) < 2) lvl <- groups_in_data
coldata$group <- factor(coldata$group, levels = lvl)

if (nlevels(coldata$group) < 2) {
  stop("group 只有一个水平：", as.character(unique(coldata$group)), "；无法做 PLL vs OPLL 对比。")
}

dds <- DESeqDataSetFromMatrix(
  countData = round(pb_counts),
  colData   = as.data.frame(coldata) %>% tibble::column_to_rownames(".agg_id"),
  design    = ~ group
)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds, quiet = TRUE)

lvl_now <- levels(coldata$group)
contrast_vec <- c("group", lvl_now[min(2, length(lvl_now))], lvl_now[1])

res <- results(dds, contrast = contrast_vec)
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(padj)

write.csv(res_df, file = "DESeq2_PLLvsOPLL_allGenes.csv", row.names = FALSE)

res_genes <- res_df %>% filter(gene %in% genes_present)
write.csv(res_genes, file = "DESeq2_PLLvsOPLL_TGFB1_IL1B.csv", row.names = FALSE)

norm_counts <- counts(dds, normalized = TRUE) %>% as.matrix()
norm_counts <- norm_counts[rownames(norm_counts) %in% genes_present, , drop = FALSE]

plot_df <- norm_counts %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = ".agg_id", values_to = "norm_count") %>%
  left_join(coldata, by = ".agg_id") %>%
  mutate(
    log2_norm = log2(norm_count + 1),
    gene = factor(gene, levels = genes_present),
    group = factor(group, levels = levels(coldata$group)),
    x_num = as.numeric(group)
  )

write.csv(plot_df, file = "pseudobulk_normalizedCounts_TGFB1_IL1B_long.csv", row.names = FALSE)

p_to_star <- function(p, cuts = c(0.001, 0.01, 0.05)) {
  if (is.na(p)) return("NA")
  if (p < cuts[1]) return("***")
  if (p < cuts[2]) return("**")
  if (p < cuts[3]) return("*")
  return("ns")
}

sig_df <- res_genes %>%
  mutate(
    p_use = ifelse(is.na(padj), pvalue, padj),
    star  = vapply(p_use, p_to_star, character(1), cuts = SIG_CUTS)
  ) %>%
  select(gene, star)

y_df <- plot_df %>%
  group_by(gene, celltype) %>%
  summarise(ymax = max(log2_norm, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    y_bracket = ymax + 0.25,
    y_text    = ymax + 0.35
  )

ann_df <- y_df %>%
  left_join(sig_df, by = "gene") %>%
  mutate(
    x1 = 1,
    x2 = 2,
    tick = 0.06
  )

if (length(levels(plot_df$group)) != 2) {
  warning("group 水平不是2个，跳过显著性标注。")
}

p1 <- ggplot(plot_df, aes(x = x_num, y = log2_norm)) +
  geom_boxplot(aes(group = group), outlier.shape = NA, width = 0.55) +
  geom_jitter(aes(group = group), width = 0.10, height = 0, size = 2, alpha = 0.80) +
  facet_grid(gene ~ celltype, scales = "free_y") +
  scale_x_continuous(breaks = c(1, 2), labels = levels(plot_df$group)) +
  labs(
    title = paste0("Pseudobulk expression: TGFB1 & IL1B (", ASSAY_USE, "/", SLOT_USE, ")"),
    subtitle = paste0("DESeq2 normalized counts; contrast: ", contrast_vec[2], " vs ", contrast_vec[3]),
    x = NULL, y = "log2(normalized counts + 1)"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

if (length(levels(plot_df$group)) == 2) {
  p1 <- p1 +
    geom_segment(
      data = ann_df,
      aes(x = x1, xend = x2, y = y_bracket, yend = y_bracket),
      inherit.aes = FALSE,
      linewidth = 0.6
    ) +
    geom_segment(
      data = ann_df,
      aes(x = x1, xend = x1, y = y_bracket - tick, yend = y_bracket),
      inherit.aes = FALSE,
      linewidth = 0.6
    ) +
    geom_segment(
      data = ann_df,
      aes(x = x2, xend = x2, y = y_bracket - tick, yend = y_bracket),
      inherit.aes = FALSE,
      linewidth = 0.6
    ) +
    geom_text(
      data = ann_df,
      aes(x = 1.5, y = y_text, label = star),
      inherit.aes = FALSE,
      size = 4,
      vjust = 0
    )
}

ggsave("pseudobulk_TGFB1_IL1B_log2norm_boxplot_sig.pdf", p1, width = 9, height = 5)
ggsave("pseudobulk_TGFB1_IL1B_log2norm_boxplot_sig.png", p1, width = 9, height = 5, dpi = 300)

cat("\n================= DONE =================\n")
cat("WORKDIR: ", normalizePath(WORKDIR, winslash = "/", mustWork = FALSE), "\n")
cat("Assay/slot: ", ASSAY_USE, "/", SLOT_USE, "\n")
cat("SAMPLE_COL: ", SAMPLE_COL, "\n")
cat("GROUP_COL : ", GROUP_COL, "\n")
cat("Celltype pseudobulk: ", DO_CELLTYPE_PSEUDOBULK, "\n")
cat("Contrast: ", contrast_vec[2], " vs ", contrast_vec[3], "\n\n")

cat("Target gene DESeq2 results (if present):\n")
print(res_genes)

cat("\nFiles:\n")
cat(" - DESeq2_PLLvsOPLL_allGenes.csv\n")
cat(" - DESeq2_PLLvsOPLL_TGFB1_IL1B.csv\n")
cat(" - pseudobulk_normalizedCounts_TGFB1_IL1B_long.csv\n")
cat(" - pseudobulk_TGFB1_IL1B_log2norm_boxplot_sig.(pdf/png)\n")
