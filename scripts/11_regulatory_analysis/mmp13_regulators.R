source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))

library(dplyr)

REG_FILE <- project_path("data_root", "auxiliary/oa/regulons/output.reg.csv")

message("正在读取 regulon 文件...")

reg_data <- read.csv(REG_FILE, stringsAsFactors = FALSE, check.names = FALSE)

rows_with_mmp13 <- apply(reg_data, 1, function(row) {
  grepl("MMP13", paste(row, collapse = " "))
})

target_rows <- reg_data[rows_with_mmp13, ]

if (nrow(target_rows) > 0) {

  message(paste0("\n>>> 成功！共发现 ", nrow(target_rows), " 个 TF 调控网络包含 MMP13。"))

  raw_tfs <- target_rows[, 1]
  clean_tfs <- gsub("\\s*\\(\\+\\)", "", raw_tfs)
  clean_tfs <- gsub("_extended", "", clean_tfs)

  result_df <- data.frame(
    Original_Name = raw_tfs,
    Clean_Name = clean_tfs,
    Regulon_Size = target_rows[, 2]
  )

  print("========== 调控 MMP13 的 TF 列表 ==========")
  print(result_df)

  if (any(grepl("ETV1", clean_tfs))) {
    message("\n★ ★ ★ 重大发现 (JACKPOT) ★ ★ ★")
    message("ETV1 被直接鉴定为 MMP13 的上游调控因子！")
    message("你的科研链条逻辑（生信部分）已经闭环。")
  } else {
    message("\n分析：ETV1 未直接出现。请检查列表中是否有 AP-1 家族 (JUN, FOS) 或 NFkB 家族。")
    message("如果出现了 JUN/FOS，说明是 ETV1-AP1 协同调控模式。")
  }

  write.csv(result_df, file.path(dirname(REG_FILE), "MMP13_Upstream_TFs.csv"), row.names = FALSE)

} else {
  message("错误：虽然刚才 grep 说有，但提取行时失败了。请检查文件编码或格式。")
}
