source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))

library(Seurat)
library(dplyr)
library(ggplot2)
library(fgsea)

seurat_object <- readRDS(project_path("data_root", "auxiliary/oa/annotated_seurat.rds"))
obj <- seurat_object
DefaultAssay(obj) <- "RNA"

map_genes_case_insensitive <- function(gene_vec, all_genes) {
  gene_vec_unique <- unique(gene_vec)
  idx        <- match(toupper(gene_vec_unique), toupper(all_genes))
  idx_valid  <- idx[!is.na(idx)]
  mapped     <- all_genes[idx_valid]
  unique(mapped)
}

safe_file_name <- function(x) {
  x <- gsub("[\\\\/:*?\"<>|]", "_", x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

score_from_gmt_dir <- function(obj,
                               dir_path,
                               prefix   = "set",
                               plot_dir = NULL,
                               width_px  = 1091,
                               height_px = 827,
                               dpi       = 300) {
  all_genes <- rownames(obj)

  if (is.null(plot_dir)) {
    plot_dir <- dir_path
  }
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }

  width_in  <- width_px  / dpi
  height_in <- height_px / dpi

  gmt_files <- list.files(
    path       = dir_path,
    pattern    = "\\.gmt$",
    full.names = TRUE
  )

  if (length(gmt_files) == 0) {
    stop("目录中没有找到 .gmt 文件: ", dir_path)
  }

  cat("在目录中找到", length(gmt_files), "个 gmt 文件：\n", dir_path, "\n\n")

  for (f in gmt_files) {
    cat(">>> 处理 GMT 文件:", basename(f), "\n")

    pathways_raw <- fgsea::gmtPathways(f)
    gs_names_all <- names(pathways_raw)
    cat("    该文件中共有基因集:", length(gs_names_all), "\n")

    pathways_mapped_all <- lapply(
      pathways_raw,
      map_genes_case_insensitive,
      all_genes = all_genes
    )
    genes_count    <- sapply(pathways_mapped_all, length)
    enough_genes   <- genes_count >= 3

    cat("    命中 ≥ 3 个基因的基因集数:", sum(enough_genes), "\n")

    score_cols_expected_all <- paste0(prefix, "_", make.names(gs_names_all))
    existing_cols           <- colnames(obj@meta.data)
    already_has_col         <- score_cols_expected_all %in% existing_cols

    if (any(already_has_col)) {
      cat("    以下基因集得分列已存在，将跳过重新计算：\n")
      print(gs_names_all[already_has_col])
      cat("\n")
    }

    if (any(!enough_genes)) {
      cat("    以下基因集命中基因 < 3，无法计算得分：\n")
      print(gs_names_all[!enough_genes])
      cat("\n")
    }

    to_score <- enough_genes & !already_has_col

    if (any(to_score)) {
      cat("    需要新增打分的基因集数:", sum(to_score), "\n")

      features_for_score <- pathways_mapped_all[to_score]
      old_cols <- colnames(obj@meta.data)

      set.seed(123)
      obj <- AddModuleScore(
        object   = obj,
        features = features_for_score,
        name     = paste0(prefix, "_tmp")
      )

      new_cols <- setdiff(colnames(obj@meta.data), old_cols)

      if (length(new_cols) != sum(to_score)) {
        warning("新列数量 (", length(new_cols),
                ") 与实际新增基因集数量 (", sum(to_score),
                ") 不一致，请检查 AddModuleScore 行为。")
      }

      n_to_rename          <- min(length(new_cols), sum(to_score))
      expected_new_names   <- score_cols_expected_all[to_score][seq_len(n_to_rename)]

      for (j in seq_len(n_to_rename)) {
        colnames(obj@meta.data)[colnames(obj@meta.data) == new_cols[j]] <- expected_new_names[j]
      }
    } else {
      cat("    没有需要新增打分的基因集。\n")
    }

    existing_cols_after <- colnames(obj@meta.data)
    plot_idx <- enough_genes & (score_cols_expected_all %in% existing_cols_after)

    cat("    准备绘图的基因集数:", sum(plot_idx), "\n\n")

    for (i in which(plot_idx)) {
      gs_name   <- gs_names_all[i]
      score_col <- score_cols_expected_all[i]

      p <- FeaturePlot(
        obj,
        features = score_col,
        pt.size  = 1,
        alpha    = 0.8,
        cols     = c("#66cbdf", "#cbdf66", "#df66cb")
      ) + ggtitle(paste0(prefix, " : ", gs_name))

      file_pdf <- file.path(
        plot_dir,
        paste0(safe_file_name(gs_name), "_featureplot.pdf")
      )
      if (file.exists(file_pdf)) next

      pdf(
        file   = file_pdf,
        width  = width_in,
        height = height_in
      )
      print(p)
      dev.off()
    }

    cat("    已完成文件:", basename(f), " 的打分和绘图。\n\n")
  }

  return(obj)
}

deg_dir <- project_path("data_root", "auxiliary/oa/gene_sets/deg")
syn_dir <- project_path("data_root", "auxiliary/oa/gene_sets/syn")

obj <- score_from_gmt_dir(
  obj,
  dir_path = deg_dir,
  prefix   = "deg",
  plot_dir = deg_dir,
  width_px  = 3452,
  height_px = 2293,
  dpi       = 300
)

obj <- score_from_gmt_dir(
  obj,
  dir_path = syn_dir,
  prefix   = "syn",
  plot_dir = syn_dir,
  width_px  = 3452,
  height_px = 2293,
  dpi       = 300
)

deg_cols <- grep("^deg_", colnames(obj@meta.data), value = TRUE)
syn_cols <- grep("^syn_", colnames(obj@meta.data), value = TRUE)

cat("当前降解相关得分列（deg_ 开头）数量:", length(deg_cols), "\n")
cat("当前合成相关得分列（syn_ 开头）数量:", length(syn_cols), "\n")

head(deg_cols)
head(syn_cols)
