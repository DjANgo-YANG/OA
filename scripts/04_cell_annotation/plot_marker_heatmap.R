rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
while (grDevices::dev.cur() > 1) grDevices::dev.off()
options(stringsAsFactors = FALSE)

WORKDIR    <- project_path("results_root", "marker_heatmaps")
SEURAT_RDS <- project_path("data_root", "auxiliary/marker_seurat.rds")

genes_use <- c(
  "COL10A1","MMP13","IBSP","SPP1","RUNX2","IHH","MEF2C","ALPL","VEGFA"
)

ASSAY_USE <- NULL
SLOT_USE  <- "data"
GROUP_BY  <- "cell_type"

OUT_PREFIX <- "Heatmap_HypertrophicMarkers_RowZscore"

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

if (!dir.exists(WORKDIR)) dir.create(WORKDIR, recursive = TRUE)
setwd(WORKDIR)

stopifnot(file.exists(SEURAT_RDS))
obj <- readRDS(SEURAT_RDS)
stopifnot(inherits(obj, "Seurat"))

if (!is.null(ASSAY_USE)) {
  if (!ASSAY_USE %in% Assays(obj)) {
    stop(sprintf("ASSAY_USE '%s' 不在对象 assays 中：%s", ASSAY_USE, paste(Assays(obj), collapse = ", ")))
  }
  DefaultAssay(obj) <- ASSAY_USE
}

meta_cols <- colnames(obj@meta.data)
group_col <- NULL
if (!is.null(GROUP_BY) && GROUP_BY %in% meta_cols) {
  group_col <- GROUP_BY
} else if ("cell_type" %in% meta_cols) {
  group_col <- "cell_type"
} else if ("seurat_clusters" %in% meta_cols) {
  group_col <- "seurat_clusters"
} else {
  obj@meta.data[[".idents_fallback"]] <- as.character(Idents(obj))
  group_col <- ".idents_fallback"
}
obj@meta.data[[group_col]] <- factor(as.character(obj@meta.data[[group_col]]))

genes_upper <- toupper(genes_use)
obj_genes_upper <- toupper(rownames(obj))
present_flag <- genes_upper %in% obj_genes_upper

genes_present <- genes_use[present_flag]
genes_missing <- genes_use[!present_flag]

if (length(genes_missing) > 0) {
  message("[WARN] 这些基因不在对象 rownames 中，将被跳过：\n", paste(genes_missing, collapse = ", "))
}
if (length(genes_present) == 0) {
  stop("所有输入基因都不在对象 rownames 中：请检查 gene symbol vs Ensembl ID / assay / 物种。")
}

avg_list <- Seurat::AverageExpression(
  object   = obj,
  features = genes_present,
  assays   = DefaultAssay(obj),
  slot     = SLOT_USE,
  group.by = group_col,
  verbose  = FALSE
)

mat <- avg_list[[ DefaultAssay(obj) ]]
mat <- mat[genes_present, , drop = FALSE]

row_z <- t(apply(mat, 1, function(x) {
  s <- stats::sd(x)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x)) / s
}))
rownames(row_z) <- rownames(mat)
colnames(row_z) <- colnames(mat)

df <- as.data.frame(row_z) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "group",
    values_to = "z"
  )

df$gene  <- factor(df$gene, levels = rev(genes_present))
df$group <- factor(df$group, levels = colnames(row_z))

p <- ggplot(df, aes(x = group, y = gene, fill = z)) +
  geom_tile(color = "white", linewidth = 0.3) +
  theme_bw() +
  theme(
    axis.title  = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid  = element_blank()
  ) +
  labs(
    title = sprintf("Row Z-score heatmap (avg %s, group.by=%s, assay=%s)",
                    SLOT_USE, group_col, DefaultAssay(obj))
  ) +
  scale_fill_gradient(
    low  = "lightyellow",
    high = "red",
    name = "Row Z-score",
    oob  = scales::squish
  )

print(p)

ggsave(
  filename = paste0(OUT_PREFIX, ".pdf"),
  plot     = p,
  width    = max(6, nlevels(df$group) * 0.6),
  height   = max(3, nlevels(df$gene)  * 0.35),
  units    = "in"
)
ggsave(
  filename = paste0(OUT_PREFIX, ".png"),
  plot     = p,
  width    = max(6, nlevels(df$group) * 0.6),
  height   = max(3, nlevels(df$gene)  * 0.35),
  units    = "in",
  dpi      = 300
)

message("[DONE] 输出：", paste0(OUT_PREFIX, ".pdf / .png"))
