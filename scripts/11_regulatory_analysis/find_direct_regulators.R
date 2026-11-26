source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))

library(dplyr)
library(stringr)

safe_write_csv <- function(data, base_path) {
  final_path <- base_path

  if (file.exists(base_path)) {
    tryCatch({

      con <- file(base_path, open = "r+")
      close(con)
    }, error = function(e) {

      timestamp <- format(Sys.time(), "%H%M%S")
      final_path <<- sub("\\.csv$", paste0("_v", timestamp, ".csv"), base_path)
      message("⚠️ 原文件被占用(可能是Excel打开了)，将自动保存为: ", basename(final_path))
    })
  }
  write.csv(data, final_path, row.names = FALSE)
  return(final_path)
}

clean_tf_name <- function(x) {
  x <- as.character(x)
  x <- gsub("\\s*\\(\\+\\)", "", x)
  x <- gsub("_extended", "", x)
  x <- gsub("^\\s+|\\s+$", "", x)
  x <- gsub(".*__", "", x)
  return(x)
}

COR_FILE <- project_path("results_root", "tf_targets/TF_Pseudotime_Correlation.csv")
REG_FILE <- project_path("data_root", "auxiliary/oa/regulons/output.reg.csv")
OUT_DIR  <- project_path("data_root", "auxiliary/oa/regulons")

if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR)

if (!file.exists(COR_FILE)) stop("❌ 找不到相关性文件: ", COR_FILE)

message(">>> 正在读取相关性数据...")
cor_df <- read.csv(COR_FILE, stringsAsFactors = FALSE)

first_col_is_numeric <- is.numeric(cor_df[[1]]) || all(grepl("^[0-9]+$", cor_df[[1]][1:5]))
tf_col_idx <- if (first_col_is_numeric) 2 else 1

raw_cor_tfs <- cor_df[[tf_col_idx]]
clean_cor_tfs <- clean_tf_name(raw_cor_tfs)

message(paste0("   - 识别到 ", length(unique(clean_cor_tfs)), " 个相关性 TF (例如: ", head(clean_cor_tfs, 1), ")"))

if (!file.exists(REG_FILE)) stop("❌ 找不到 SCENIC 结果文件: ", REG_FILE)

message(">>> 正在读取 SCENIC Regulon 数据...")

reg_lines <- readLines(REG_FILE)

message(">>> 开始全景扫描 MMP13 调控关系...")

results_list <- list()

for (i in 1:length(reg_lines)) {
  line_content <- reg_lines[i]

  if (grepl("MMP13", line_content, fixed = TRUE)) {

    parts <- strsplit(line_content, "[,;\t]")[[1]]

    candidate_tf <- parts[1]

    if (grepl("^[0-9]+$", candidate_tf)) candidate_tf <- parts[2]

    current_tf_clean <- clean_tf_name(candidate_tf)

    is_correlated <- current_tf_clean %in% clean_cor_tfs

    results_list[[length(results_list) + 1]] <- data.frame(
      TF_Raw = candidate_tf,
      TF_Clean = current_tf_clean,
      Target = "MMP13",
      Is_Pseudotime_Correlated = is_correlated,
      Raw_Line_Index = i,
      stringsAsFactors = FALSE
    )
  }
}

if (length(results_list) > 0) {
  final_result <- do.call(rbind, results_list)

  final_result <- unique(final_result)

  message("\n========================================================")
  message("✅ 分析完成！捕获到以下调控 MMP13 的转录因子：")
  message("========================================================")

  hits <- final_result[final_result$Is_Pseudotime_Correlated == TRUE, ]
  others <- final_result[final_result$Is_Pseudotime_Correlated == FALSE, ]

  if (nrow(hits) > 0) {
    message("\n[重点关注] 既调控 MMP13 又随拟时序显著变化的 TF：")
    print(hits[, c("TF_Clean", "Raw_Line_Index")])

    if ("ETV1" %in% hits$TF_Clean) {
      message("\n★ 确认：ETV1 在列表中！科研链条闭环成功。")
    }
  } else {
    message("\n[提示] 找到了调控 MMP13 的 TF，但它们似乎不在拟时序相关列表中。")
    message("请检查清洗后的 TF 名字是否匹配。")
  }

  if (nrow(others) > 0) {
    message("\n[参考] 其他调控 MMP13 的 TF (非拟时序显著)：")
    print(others$TF_Clean)
  }

  save_path <- file.path(OUT_DIR, "MMP13_Direct_Regulators_Robust.csv")
  actual_path <- safe_write_csv(final_result, save_path)
  message("\n结果已保存至: ", actual_path)

} else {
  message("\n⚠️ 警告：在 SCENIC 文件中完全未扫描到 'MMP13' 字符串。")
  message("可能性 1: SCENIC 真的没算出来 (False Negative)。")
  message("可能性 2: 基因名大小写不同 (如 mmp13)。")
}
