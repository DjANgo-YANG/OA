rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
while (grDevices::dev.cur() > 1) grDevices::dev.off()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(DESeq2)
})

MRNA1_COUNTS <- project_path("data_root", "auxiliary/trajectory/pseudobulk/GSE69787_RNA_All_Counts.txt")
MRNA2_COUNTS <- project_path("data_root", "auxiliary/trajectory/pseudobulk/GSE183752_All.counts.exp.txt")
stopifnot(file.exists(MRNA1_COUNTS), file.exists(MRNA2_COUNTS))

WORKDIR <- project_path("data_root", "auxiliary/trajectory/pseudobulk/BULK_TGFB1_IL1B_EXPLICIT_MAPPING_SIG")
dir.create(WORKDIR, showWarnings = FALSE, recursive = TRUE)
setwd(WORKDIR)

GENES <- c("TGFB1", "IL1B")

GSE69787_PLL  <- c("PLL-1", "PLL-2", "PLL-3")
GSE69787_OPLL <- c("OPLL-1", "OPLL-2", "OPLL-3")

GSE183752_OPLL <- c("OPLL-A-hTERT", "OPLL-B-SV40", "OPLL-C-SV40", "OPLL-NC")
GSE183752_PLL  <- c("PLL-A-hTERT",  "PLL-B-SV40",  "PLL-C-SV40",  "PLL-NC")

SIG_CUTS <- c(0.001, 0.01, 0.05)

read_counts_matrix <- function(path) {

  dt <- suppressWarnings(data.table::fread(path, data.table = FALSE, check.names = FALSE))
  if (ncol(dt) < 2) stop("文件列数<2，无法识别counts矩阵：", path)

  feature_col <- colnames(dt)[1]
  features <- dt[[feature_col]]
  dt[[feature_col]] <- NULL

  is_num_col <- vapply(dt, function(x) {
    suppressWarnings({
      xx <- as.numeric(x)
      sum(!is.na(xx)) > 0
    })
  }, logical(1))
  if (!any(is_num_col)) stop("未检测到可用的数值列（counts列）：", path)

  dt_num <- dt[, is_num_col, drop = FALSE]
  mat <- as.matrix(sapply(dt_num, function(x) suppressWarnings(as.numeric(x))))
  rownames(mat) <- make.unique(as.character(features))
  colnames(mat) <- colnames(dt_num)

  mat[is.na(mat)] <- 0
  mat <- round(mat)
  storage.mode(mat) <- "integer"
  mat
}

validate_and_build_coldata <- function(samples_in_mat, pll_vec, opll_vec, dataset_key, outdir) {

  expected <- c(pll_vec, opll_vec)

  missing <- setdiff(expected, samples_in_mat)
  extra   <- setdiff(samples_in_mat, expected)

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  writeLines(paste(samples_in_mat, collapse = "\n"), con = file.path(outdir, "ACTUAL_COLNAMES.txt"))
  writeLines(paste(expected, collapse = "\n"),       con = file.path(outdir, "EXPECTED_COLNAMES.txt"))

  if (length(missing) > 0 || length(extra) > 0) {
    stop(
      dataset_key, "：列名与显式 mapping 不一致。\n",
      "缺失列名(missing): ", ifelse(length(missing)>0, paste(missing, collapse = ", "), "无"), "\n",
      "多余列名(extra):   ", ifelse(length(extra)>0,   paste(extra, collapse = ", "),   "无"), "\n",
      "请查看：\n",
      " - ", file.path(outdir, "ACTUAL_COLNAMES.txt"), "\n",
      " - ", file.path(outdir, "EXPECTED_COLNAMES.txt"), "\n"
    )
  }

  coldata <- data.frame(
    sample = expected,
    group  = c(rep("PLL", length(pll_vec)), rep("OPLL", length(opll_vec))),
    stringsAsFactors = FALSE
  )

  coldata$group <- factor(coldata$group, levels = c("PLL", "OPLL"))
  coldata
}

p_to_star <- function(p, cuts = c(0.001, 0.01, 0.05)) {

  if (is.na(p)) return("NA")
  if (p < cuts[1]) return("***")
  if (p < cuts[2]) return("**")
  if (p < cuts[3]) return("*")
  return("ns")
}

fmt_p <- function(p) {

  if (is.na(p)) return("NA")
  if (p < 1e-4) return(format(p, scientific = TRUE, digits = 2))
  return(sprintf("%.4f", p))
}

run_bulk_explicit_mapping <- function(counts_mat, dataset_key, outdir, genes_of_interest, pll_vec, opll_vec) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  samples <- colnames(counts_mat)
  coldata <- validate_and_build_coldata(
    samples_in_mat = samples,
    pll_vec = pll_vec,
    opll_vec = opll_vec,
    dataset_key = dataset_key,
    outdir = outdir
  )
  write.csv(coldata, file = file.path(outdir, "sample_coldata_USED.csv"), row.names = FALSE)

  counts_mat <- counts_mat[, coldata$sample, drop = FALSE]

  genes_present <- intersect(genes_of_interest, rownames(counts_mat))
  genes_missing <- setdiff(genes_of_interest, genes_present)

  if (length(genes_missing) > 0) {
    writeLines(
      paste0("以下目标基因在该矩阵中未找到：", paste(genes_missing, collapse = ", ")),
      con = file.path(outdir, "TARGET_GENES_MISSING.txt")
    )
  }
  if (length(genes_present) > 0) {
    raw_df <- counts_mat[genes_present, , drop = FALSE] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("gene")
    write.csv(raw_df, file = file.path(outdir, "TGFB1_IL1B_rawCounts.csv"), row.names = FALSE)
  }

  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData   = tibble::column_to_rownames(coldata, "sample"),
    design    = ~ group
  )
  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds, quiet = TRUE)

  contrast_vec <- c("group", "OPLL", "PLL")

  res <- results(dds, contrast = contrast_vec)
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(padj)

  write.csv(res_df, file = file.path(outdir, "DESeq2_allFeatures.csv"), row.names = FALSE)

  res_goi <- res_df %>% filter(gene %in% genes_present)
  if (nrow(res_goi) > 0) {
    write.csv(res_goi, file = file.path(outdir, "DESeq2_TGFB1_IL1B.csv"), row.names = FALSE)
  }

  if (length(genes_present) > 0) {

    norm_counts <- counts(dds, normalized = TRUE) %>% as.matrix()
    norm_counts <- norm_counts[genes_present, , drop = FALSE]

    plot_df <- norm_counts %>%
      as.data.frame() %>%
      tibble::rownames_to_column("gene") %>%
      pivot_longer(cols = -gene, names_to = "sample", values_to = "norm_count") %>%
      left_join(coldata, by = "sample") %>%
      mutate(
        log2_norm = log2(norm_count + 1),
        gene = factor(gene, levels = genes_present),
        group = factor(group, levels = c("PLL", "OPLL")),

        x_num = ifelse(group == "PLL", 1, 2)
      )

    write.csv(plot_df, file = file.path(outdir, "TGFB1_IL1B_normalizedCounts_long.csv"), row.names = FALSE)

    sig_df <- res_goi %>%
      mutate(
        p_use = ifelse(is.na(padj), pvalue, padj),
        star  = vapply(p_use, p_to_star, character(1), cuts = SIG_CUTS),
        ptxt  = vapply(p_use, fmt_p, character(1)),
        label = paste0(star, "\n", "padj=", ptxt)
      ) %>%
      select(gene, p_use, label)

    y_df <- plot_df %>%
      group_by(gene) %>%
      summarise(ymax = max(log2_norm, na.rm = TRUE), .groups = "drop") %>%
      mutate(
        y_bracket = ymax + 0.25,
        y_text    = ymax + 0.35
      )

    ann_df <- y_df %>%
      left_join(sig_df, by = "gene") %>%
      mutate(

        x1 = 1, x2 = 2,

        tick = 0.06
      )

    p <- ggplot(plot_df, aes(x = x_num, y = log2_norm)) +
      geom_boxplot(aes(group = group), outlier.shape = NA, width = 0.55) +
      geom_jitter(aes(group = group), width = 0.10, height = 0, size = 2, alpha = 0.85) +
      facet_wrap(~ gene, scales = "free_y", ncol = 2) +
      scale_x_continuous(breaks = c(1, 2), labels = c("PLL", "OPLL")) +
      labs(
        title = paste0(dataset_key, ": TGFB1 & IL1B (bulk mRNA)"),
        subtitle = "DESeq2 normalized counts; contrast: OPLL vs PLL",
        x = NULL, y = "log2(normalized counts + 1)"
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.minor = element_blank()
      )

    p <- p +
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
        aes(x = 1.5, y = y_text, label = label),
        inherit.aes = FALSE,
        size = 3.5,
        lineheight = 0.95
      )

    ggsave(file.path(outdir, "TGFB1_IL1B_log2norm_boxplot_sig.pdf"), p, width = 9, height = 5)
    ggsave(file.path(outdir, "TGFB1_IL1B_log2norm_boxplot_sig.png"), p, width = 9, height = 5, dpi = 300)
  }

  invisible(NULL)
}

mat1 <- read_counts_matrix(MRNA1_COUNTS)
run_bulk_explicit_mapping(
  counts_mat = mat1,
  dataset_key = "GSE69787",
  outdir = file.path(WORKDIR, "GSE69787"),
  genes_of_interest = GENES,
  pll_vec  = GSE69787_PLL,
  opll_vec = GSE69787_OPLL
)

mat2 <- read_counts_matrix(MRNA2_COUNTS)
run_bulk_explicit_mapping(
  counts_mat = mat2,
  dataset_key = "GSE183752",
  outdir = file.path(WORKDIR, "GSE183752"),
  genes_of_interest = GENES,
  pll_vec  = GSE183752_PLL,
  opll_vec = GSE183752_OPLL
)

cat("\n================= DONE =================\n")
cat("WORKDIR:\n  ", normalizePath(WORKDIR, winslash = "/", mustWork = FALSE), "\n\n")
cat("输出目录:\n")
cat(" - ", normalizePath(file.path(WORKDIR, "GSE69787"), winslash = "/", mustWork = FALSE), "\n")
cat(" - ", normalizePath(file.path(WORKDIR, "GSE183752"), winslash = "/", mustWork = FALSE), "\n")
cat("\n图文件（带显著性）示例：TGFB1_IL1B_log2norm_boxplot_sig.pdf/png\n")
