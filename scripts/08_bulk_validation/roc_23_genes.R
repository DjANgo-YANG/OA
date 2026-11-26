rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pROC)
})

sample_info_file <- project_path("data_root", "bulk/sample_info.csv")
tpm_file         <- project_path("data_root", "bulk/tpm.csv")
outdir <- Sys.getenv("OAGLY_ROC_OUTDIR", unset = project_path("data_root", "bulk/AUC_23genes"))

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(sample_info_file))
stopifnot(file.exists(tpm_file))

genes_23 <- c(
  "MDH1","DCN","PKM","LDHA","SLC35A3","PLOD2","CD44","ENO1","PPIA","AK3",
  "TGFBI","AGL","NASP","LHX9","COPB2","PAM","VCAN","ELF3","SDC2","ME2",
  "PGK1","HDLBP","COG2"
)

sample_info <- fread(sample_info_file)

cat("sample_info columns:\n")
print(colnames(sample_info))

sample_col_candidates <- c("sample", "sample_id", "Sample", "SampleID", "sampleID")
sample_col <- sample_col_candidates[sample_col_candidates %in% colnames(sample_info)][1]

if (is.na(sample_col)) {
  stop("sample_info.csv 中未找到样本列。请检查是否包含以下列之一: ",
       paste(sample_col_candidates, collapse = ", "))
}

if (!"condition" %in% colnames(sample_info)) {
  stop("sample_info.csv 中不存在列: condition")
}

sample_info <- sample_info %>%
  dplyr::select(sample = all_of(sample_col), condition)

sample_info <- sample_info %>%
  mutate(
    condition = as.character(condition),
    group = case_when(
      toupper(condition) == "HC" ~ "HC",
      toupper(condition) == "OA" ~ "OA",
      toupper(condition) == "CONTROL" ~ "HC",
      toupper(condition) %in% c("WB", "NWB") ~ "OA",
      grepl("^HC", condition, ignore.case = TRUE) ~ "HC",
      grepl("^OA", condition, ignore.case = TRUE) ~ "OA",
      TRUE ~ NA_character_
    )
  )

cat("\nRaw condition table:\n")
print(table(sample_info$condition, useNA = "ifany"))

cat("\nConverted group table:\n")
print(table(sample_info$group, useNA = "ifany"))

sample_info <- sample_info %>% filter(group %in% c("HC", "OA"))

if (nrow(sample_info) == 0) {
  stop("condition 列未能成功转换出 HC / OA，请检查 sample_info$condition 的具体内容。")
}

tpm <- fread(tpm_file)

cat("\ntpm columns (head):\n")
print(head(colnames(tpm), 10))

colnames(tpm)[1] <- "gene"

tpm <- tpm %>%
  distinct(gene, .keep_all = TRUE)

tpm_sub <- tpm %>% filter(gene %in% genes_23)

genes_present <- intersect(genes_23, tpm_sub$gene)
genes_missing <- setdiff(genes_23, tpm_sub$gene)

writeLines(genes_present, file.path(outdir, "genes_present.txt"))
writeLines(genes_missing, file.path(outdir, "genes_missing.txt"))

cat("\nGenes present:", length(genes_present), "\n")
if (length(genes_missing) > 0) {
  cat("Genes missing:\n")
  print(genes_missing)
}

if (length(genes_missing) > 0) {
  stop("The TPM matrix lacks candidate genes: ", paste(genes_missing, collapse = ", "))
}

expr_sample_cols <- setdiff(colnames(tpm_sub), "gene")
common_samples <- intersect(expr_sample_cols, sample_info$sample)

cat("\nMatched samples:", length(common_samples), "\n")
print(common_samples)

if (length(common_samples) < 4) {
  stop("表达矩阵和 sample_info 匹配到的样本太少。")
}

sample_info2 <- sample_info %>%
  filter(sample %in% common_samples) %>%
  arrange(match(sample, common_samples))

expr_mat <- as.data.frame(tpm_sub)
rownames(expr_mat) <- expr_mat$gene
expr_mat$gene <- NULL
expr_mat <- expr_mat[, sample_info2$sample, drop = FALSE]

expr_log <- log2(as.matrix(expr_mat) + 1)

group_factor <- factor(sample_info2$group, levels = c("HC", "OA"))

auc_res <- lapply(rownames(expr_log), function(g) {
  x <- as.numeric(expr_log[g, ])

  roc_obj <- roc(
    response = group_factor,
    predictor = x,
    levels = c("HC", "OA"),
    direction = "<",
    quiet = TRUE
  )

  ci_auc <- ci.auc(roc_obj)

  data.frame(
    Gene = g,
    AUC = as.numeric(auc(roc_obj)),
    CI_low = as.numeric(ci_auc[1]),
    CI_high = as.numeric(ci_auc[3]),
    stringsAsFactors = FALSE
  )
})

auc_tbl <- bind_rows(auc_res) %>%
  arrange(desc(AUC))

write.csv(auc_tbl, file.path(outdir, "AUC_23genes_results.csv"), row.names = FALSE)

cat("\nAUC results:\n")
print(auc_tbl)

auc_tbl$Gene <- factor(auc_tbl$Gene, levels = auc_tbl$Gene)

p_auc <- ggplot(auc_tbl, aes(x = Gene, y = AUC)) +
  geom_col(fill = "#4e79a7", width = 0.7) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "red") +
  theme_classic(base_size = 14) +
  labs(
    title = "AUC of 22 candidate genes for distinguishing OA from HC",
    x = "Gene",
    y = "AUC"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(file.path(outdir, "AUC_23genes_barplot.pdf"), p_auc, width = 10, height = 6)
ggsave(file.path(outdir, "AUC_23genes_barplot.png"), p_auc, width = 10, height = 6, dpi = 300)

roc_dir <- file.path(outdir, "ROC_curves")
dir.create(roc_dir, showWarnings = FALSE, recursive = TRUE)

for (g in rownames(expr_log)) {
  x <- as.numeric(expr_log[g, ])

  roc_obj <- roc(
    response = group_factor,
    predictor = x,
    levels = c("HC", "OA"),
    direction = "<",
    quiet = TRUE
  )

  auc_val <- as.numeric(auc(roc_obj))

  pdf(file.path(roc_dir, paste0(g, "_ROC.pdf")), width = 6, height = 6)
  plot(roc_obj,
       col = "#1f77b4",
       lwd = 2,
       main = paste0(g, " ROC curve\nAUC = ", round(auc_val, 3)))
  abline(a = 0, b = 1, lty = 2, col = "grey50")
  dev.off()
}

expr_df <- as.data.frame(expr_log)
expr_df$Gene <- rownames(expr_df)

expr_long <- expr_df %>%
  pivot_longer(cols = -Gene, names_to = "sample", values_to = "log2TPM") %>%
  left_join(sample_info2, by = "sample")

write.csv(expr_long, file.path(outdir, "AUC_23genes_expression_long.csv"), row.names = FALSE)

p_box <- ggplot(expr_long, aes(x = group, y = log2TPM, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.8) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  theme_classic(base_size = 13) +
  labs(
    title = "Expression of 22 candidate genes in HC and OA",
    x = "Group",
    y = "log2(TPM + 1)"
  ) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  )

ggsave(file.path(outdir, "AUC_23genes_expression_boxplot.pdf"), p_box, width = 13, height = 10)
ggsave(file.path(outdir, "AUC_23genes_expression_boxplot.png"), p_box, width = 13, height = 10, dpi = 300)

top10_tbl <- auc_tbl %>% slice_head(n = 10)
write.csv(top10_tbl, file.path(outdir, "AUC_top10_genes.csv"), row.names = FALSE)

message("Done.")
message("Output dir: ", outdir)
message("Present genes: ", paste(genes_present, collapse = ", "))
if (length(genes_missing) > 0) {
  message("Missing genes: ", paste(genes_missing, collapse = ", "))
}
