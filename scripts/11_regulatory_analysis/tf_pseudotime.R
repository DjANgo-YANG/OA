rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
options(stringsAsFactors = FALSE)

SEURAT_RDS  <- project_path("data_root", "auxiliary/oa/slingshot/rds/OA_Seurat_with_slingshot.rds")

AUC_FILE    <- project_path("data_root", "auxiliary/oa/regulons/00_tables/Top231_TF_selected_for_DotPlot_by_meanAUC_allcells.csv")

OUT_DIR     <- project_path("results_root", "tf_targets")

TIME_COL    <- "pseudotime_L5"

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(data.table)

message("正在读取 Seurat 对象...")
seu <- readRDS(SEURAT_RDS)

if (!TIME_COL %in% colnames(seu@meta.data)) {
  stop(paste0("错误：在 meta.data 中找不到 ", TIME_COL, "，请检查拼写。"))
}

pt_vec <- seu@meta.data[[TIME_COL]]
names(pt_vec) <- rownames(seu@meta.data)
pt_vec <- pt_vec[!is.na(pt_vec)]

message(paste0("Lineage 5 上共有 ", length(pt_vec), " 个细胞。"))

AUC_LOOM_FILE <- project_path("data_root", "auxiliary/oa/regulons/output.auc_mtx.loom")

if (!requireNamespace("SCopeLoomR", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("aertslab/SCopeLoomR")
}

if (!requireNamespace("AUCell", quietly = TRUE)) BiocManager::install("AUCell")

library(SCopeLoomR)
library(AUCell)

message("正在读取 Loom 文件: ", AUC_LOOM_FILE)
loom <- open_loom(AUC_LOOM_FILE)

auc_data <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")

close_loom(loom)

message("读取对象的类型: ", class(auc_data))

auc_mat <- NULL

if (inherits(auc_data, "aucellResults")) {
  message("检测到 AUCell 对象，正在提取矩阵...")
  auc_mat <- AUCell::getAUC(auc_data)
} else if (is.list(auc_data) && "AUC" %in% names(auc_data)) {

  auc_mat <- auc_data$AUC
} else {

  auc_mat <- auc_data
}

auc_mat <- as.matrix(auc_mat)

message("原始矩阵维度: ", paste(dim(auc_mat), collapse = " x "))

if (nrow(auc_mat) < ncol(auc_mat)) {

  message("正在转置矩阵为 (Cells x Regulons) 格式...")
  auc_mat <- t(auc_mat)
} else {
  message("矩阵格式看似已经是 (Cells x Regulons)，跳过转置。")
}

message("最终矩阵维度: ", paste(dim(auc_mat), collapse = " x "))
message("前 3 个行名 (Cells?): ", paste(head(rownames(auc_mat), 3), collapse = ", "))
message("前 3 个列名 (Regulons?): ", paste(head(colnames(auc_mat), 3), collapse = ", "))

common_cells <- intersect(names(pt_vec), rownames(auc_mat))

if (length(common_cells) < 100) stop("错误：有拟时序和 AUC 数据的重叠细胞太少 (<100)，请检查细胞名是否一致。")

pt_use <- pt_vec[common_cells]
auc_use <- auc_mat[common_cells, ]

message("正在计算相关性 (Correlation)...")

cor_res <- apply(auc_use, 2, function(x) cor(x, pt_use, method = "spearman"))

cor_df <- data.frame(TF = names(cor_res), Correlation = cor_res)
cor_df <- cor_df %>% arrange(desc(Correlation))

write.csv(cor_df, file.path(OUT_DIR, "TF_Pseudotime_Correlation.csv"), row.names = FALSE)

top10_tfs <- head(cor_df$TF, 10)
message("与拟时序正相关性最高的 Top 10 TF:")
print(top10_tfs)

if ("ZEB1(+)" %in% top10_tfs || "ZEB1" %in% top10_tfs) message("提示: ZEB1 在 Top 10 中！")

library(grid)

p1 <- ggplot(plot_data_long, aes(x = Pseudotime, y = AUC, color = Regulon)) +
  geom_smooth(method = "loess", span = 0.7, se = TRUE, alpha = 0.1, linewidth = 1.2) +
  theme_classic() +
  scale_color_brewer(palette = "Paired") +

  labs(title = paste0("Top 10 TFs Positively Correlated with ", TIME_COL),
       x = "Pseudotime (MMP13-  --->  MMP13+)",
       y = "Regulon AUC Activity") +

  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),

    legend.position = "bottom",
    legend.direction = "horizontal",

    legend.spacing.x = unit(0.8, 'cm'),

    legend.text = element_text(size = 9),
    legend.title = element_text(face = "bold"),

    legend.margin = margin(t = 10, b = 10)
  ) +

  guides(color = guide_legend(

    ncol = 4,

    byrow = TRUE,

    keywidth = unit(1.5, "cm"),
    keyheight = unit(0.5, "cm")
  ))

ggsave(file.path(OUT_DIR, "Top10_TF_Pseudotime_Trend_Spaced.png"), p1, width = 12, height = 7, dpi = 300)
ggsave(file.path(OUT_DIR, "Top10_TF_Pseudotime_Trend_Spaced.pdf"), p1, width = 12, height = 7)

message("已生成加大间距版图表：Top10_TF_Pseudotime_Trend_Spaced.png")
