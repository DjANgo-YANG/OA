source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))

library(Seurat)
library(dplyr)

OUT_DIR <- project_path("data_root", "auxiliary/oa/regulons")
if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR)

SEURAT_FILE <- project_path("data_root", "auxiliary/oa/slingshot/rds/OA_Seurat_with_slingshot.rds")
message("正在读取 Seurat 对象: ", SEURAT_FILE)
seu <- readRDS(SEURAT_FILE)

COR_FILE <- project_path("results_root", "tf_targets/TF_Pseudotime_Correlation.csv")
message("正在读取 TF 列表: ", COR_FILE)

if(file.exists(COR_FILE)) {
  tf_data <- read.csv(COR_FILE, stringsAsFactors = FALSE)

  sig_tfs <- tf_data[, 1]
  message(">>> 成功加载 ", length(sig_tfs), " 个候选 TF")
} else {
  stop("找不到 TF 列表文件，请检查路径！")
}

mmp13_expr <- FetchData(seu, vars = "MMP13", layer = "data")
if(ncol(mmp13_expr) == 0) {

  mmp13_expr <- FetchData(seu, vars = "MMP13")
}
mmp13_vec <- mmp13_expr[, 1]

clean_tfs <- gsub("\\(\\+\\)|_extended", "", sig_tfs)

valid_tfs <- intersect(clean_tfs, rownames(seu))

message(">>> 在表达矩阵中实际找到 ", length(valid_tfs), " 个匹配的 TF")

if(length(valid_tfs) > 0) {

  tf_expr_data <- FetchData(seu, vars = valid_tfs, layer = "data")

  cor_results <- numeric(length(valid_tfs))
  names(cor_results) <- valid_tfs

  message("正在计算相关性 (Spearman)...")

  for(tf in valid_tfs) {

    tf_vec <- tf_expr_data[[tf]]

    res <- cor(tf_vec, mmp13_vec, method = "spearman")
    cor_results[tf] <- res
  }

  top_cor_tfs <- sort(cor_results, decreasing = TRUE)

  message("\n=============================================")
  message("与 MMP13 正相关性最高的 Top 10 TF:")
  message("=============================================")
  print(head(top_cor_tfs, 10))

  save_path <- file.path(OUT_DIR, "TF_MMP13_Correlation_Result.csv")
  write.csv(data.frame(TF=names(top_cor_tfs), Cor=top_cor_tfs), save_path)
  message("\n完整结果已保存至: ", save_path)

} else {
  message("错误：TF 名字匹配失败，请检查 SCENIC 输出的名字格式与 Seurat 对象中的基因名是否一致。")
}

library(ggplot2)
library(dplyr)

csv_file <- file.path(OUT_DIR, "TF_MMP13_Correlation_Result.csv")
plot_data <- read.csv(csv_file, row.names = 1)

colnames(plot_data) <- c("TF", "Cor")

top_plot_data <- plot_data %>%
  arrange(desc(Cor)) %>%
  slice_head(n = 15)

top_plot_data$Highlight <- ifelse(top_plot_data$TF == "ETV1", "Target", "Others")

p <- ggplot(top_plot_data, aes(x = reorder(TF, Cor), y = Cor, fill = Highlight)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Target" = "#E41A1C", "Others" = "#377EB8")) +
  labs(
    title = "Top TFs Correlated with MMP13 in FC",
    subtitle = "Spearman Correlation Analysis",
    x = "Transcription Factors",
    y = "Correlation Coefficient (Spearman)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(face = "bold", color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray50")
  ) +
  geom_text(aes(label = round(Cor, 2)), hjust = -0.1, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

print(p)

ggsave(filename = file.path(OUT_DIR, "Barplot_MMP13_Correlated_TFs.pdf"),
       plot = p, width = 6, height = 5)

message("绘图完成！图片已保存至: ", file.path(OUT_DIR, "Barplot_MMP13_Correlated_TFs.pdf"))
