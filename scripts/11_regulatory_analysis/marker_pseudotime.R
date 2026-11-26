rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
options(stringsAsFactors = FALSE)

GAM_RDS <- project_path("data_root", "auxiliary/oa/slingshot/rds/OA_tradeSeq_fitGAM.rds")

SEURAT_RDS <- project_path("data_root", "auxiliary/oa/slingshot/rds/OA_Seurat_with_slingshot.rds")

OUT_DIR <- project_path("results_root", "tradeseq_genes")

TARGET_GENES <- c("MMP13", "POSTN", "COL2A1")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

library(tradeSeq)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

message("正在读取 fitGAM 模型结果 (可能文件较大)...")
sce_gam <- readRDS(GAM_RDS)

message("正在读取 Seurat 对象以获取 Counts...")
seu <- readRDS(SEURAT_RDS)

if ("counts" %in% names(seu@assays$RNA)) {
  counts_mat <- GetAssayData(seu, assay = "RNA", slot = "counts")
} else {

  warning("未找到 counts slot，尝试使用 data slot...")
  counts_mat <- GetAssayData(seu, assay = "RNA", slot = "data")
}

valid_genes <- intersect(TARGET_GENES, rownames(counts_mat))
if (length(valid_genes) == 0) stop("错误：目标基因都不在表达矩阵中！请检查拼写。")

message("开始绘制基因趋势图...")

plot_list <- list()

for (gene in valid_genes) {

  message(paste0("正在绘制: ", gene))

  p <- plotSmoothers(sce_gam, counts = counts_mat, gene = gene, alpha = 0.3, border = TRUE) +
    theme_classic() +
    labs(title = gene,
         x = "Pseudotime (MMP13-  --->  MMP13+)",
         y = "Expression (Log Scale)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right") +
    scale_color_brewer(palette = "Set1", name = "Lineage")

  ggsave(file.path(OUT_DIR, paste0("Trend_", gene, ".png")), p, width = 6, height = 5)
  ggsave(file.path(OUT_DIR, paste0("Trend_", gene, ".pdf")), p, width = 6, height = 5)

  plot_list[[gene]] <- p
}

message("正在生成组合图...")

if (length(plot_list) >= 2) {
  p_combined <- wrap_plots(plot_list, ncol = 3) +
    plot_annotation(title = "Key Gene Dynamics along Pseudotime",
                    theme = theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold")))

  ggsave(file.path(OUT_DIR, "Combined_Genes_Trend.png"), p_combined, width = 15, height = 5)
  ggsave(file.path(OUT_DIR, "Combined_Genes_Trend.pdf"), p_combined, width = 15, height = 5)
}

message("========================================")
message("绘图完成！结果已保存至: ", OUT_DIR)
message("请检查 Combined_Genes_Trend.png")
