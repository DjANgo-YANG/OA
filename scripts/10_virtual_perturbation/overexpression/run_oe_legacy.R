rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
options(stringsAsFactors = FALSE)

stopifnot(getRversion() >= "4.5")

in_rds   <- project_path("results_root", "glycolysis_scoring/seurat_with_LeadingEdge_scores.rds")
out_root <- project_path("results_root", "virtual_oe_legacy")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

target_celltypes <- c("preHTC", "HTC")

min_cells <- 25
n_hvg     <- 10000
n_cells_use_cap <- 800
nCores_use <- max(1, parallel::detectCores() - 1)

force_rerun <- FALSE

gene_map <- c(
  "01_AGL" = "AGL", "02_NASP" = "NASP", "03_LHX9" = "LHX9", "04_COPB2" = "COPB2",
  "05_PAM" = "PAM", "06_VCAN" = "VCAN", "07_ELF3" = "ELF3", "08_SDC2" = "SDC2",
  "09_ME2" = "ME2", "10_PGK1" = "PGK1", "11_HDLBP" = "HDLBP", "12_MDH1" = "MDH1",
  "13_DCN" = "DCN", "14_PKM" = "PKM", "15_LDHA" = "LDHA", "16_SLC35A3" = "SLC35A3",
  "17_PLOD2" = "PLOD2", "18_CD44" = "CD44", "19_ENO1" = "ENO1", "20_PPIA" = "PPIA",
  "21_AK3" = "AK3", "22_TGFBI" = "TGFBI", "23_COG2" = "COG2"
)

cli_args <- commandArgs(trailingOnly = TRUE)
if (length(cli_args) > 0) {
  sel <- toupper(unname(gene_map)) %in% toupper(cli_args)
  if (!any(sel)) stop("命令行基因不在 gene_map 中: ", paste(cli_args, collapse = ", "))
  gene_map <- gene_map[sel]
}

log_file <- file.path(out_root, "虚拟过表达_run.log")
log_msg <- function(...) {
  msg <- sprintf("[%s] %s", Sys.time(), paste0(...))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

log_msg("=== 虚拟过表达批量运行开始 ===")
log_msg("R ", R.version.string, " | scTenifoldPerturb ",
        as.character(packageVersion("scTenifoldPerturb")))
log_msg("本次基因数: ", length(gene_map), " (", paste(unname(gene_map), collapse = ","), ")")

pkgs <- c("Matrix", "scTenifoldPerturb", "ggplot2", "dplyr", "ggrepel", "igraph", "tibble", "stringr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Missing package: ", p)
}
suppressPackageStartupMessages(lapply(pkgs, library, character.only = TRUE))
if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("Missing package: SeuratObject")

pick_meta_col <- function(obj, candidates) {
  md <- obj@meta.data
  cn <- colnames(md)
  hit <- intersect(candidates, cn)
  if (length(hit) > 0) return(hit[1])
  cn_low <- tolower(cn)
  cand_low <- tolower(candidates)
  idx <- match(cand_low, cn_low)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0) return(cn[idx[1]])
  return(NA_character_)
}

get_counts <- function(x) {
  if (inherits(x, "dgCMatrix") || is.matrix(x)) return(x)
  if (inherits(x, "Seurat")) {
    tryCatch({
      return(SeuratObject::GetAssayData(x, assay = "RNA", slot = "counts"))
    }, error = function(e) {
      return(x@assays$RNA@counts)
    })
  }
  stop("Unsupported object type for counts.")
}

fix_gene_symbol <- function(g, all_genes) {
  if (g %in% all_genes) return(g)
  tries <- unique(c(g, toupper(g), tolower(g), stringr::str_to_title(g)))
  hit <- tries[tries %in% all_genes]
  if (length(hit) > 0) return(hit[1])
  if (toupper(g) %in% all_genes) return(toupper(g))
  if (paste0(toupper(substr(g, 1, 1)), tolower(substr(g, 2, nchar(g)))) %in% all_genes) {
    return(paste0(toupper(substr(g, 1, 1)), tolower(substr(g, 2, nchar(g)))))
  }
  return(NA_character_)
}

log_msg(">>> 读取 Seurat 对象: ", in_rds)
obj <- readRDS(in_rds)
if (!inherits(obj, "Seurat")) stop("Input RDS is not a Seurat object.")

ct_col <- pick_meta_col(obj, c("cell_type", "celltype", "CellType", "cellType", "annot", "annotation", "seurat_clusters"))
if (is.na(ct_col)) stop("Cannot find cell type column in obj@meta.data. Please check metadata.")
log_msg(">>> 使用细胞类型列: ", ct_col)

md <- obj@meta.data
keep_cells <- rownames(md)[as.character(md[[ct_col]]) %in% target_celltypes]
if (length(keep_cells) < 50) {
  stop("Too few cells after filtering preHTC/HTC (n=", length(keep_cells), "). Check your cell type labels.")
}
obj_sub <- subset(obj, cells = keep_cells)

countMatrix <- get_counts(obj_sub)

missing_genes <- setdiff(unname(gene_map), rownames(countMatrix))
if (length(missing_genes) > 0) {
  stop("Missing genes in matrix: ", paste(missing_genes, collapse = ", "))
}

set.seed(123)

detected <- Matrix::rowSums(countMatrix != 0)
count_f  <- countMatrix[detected >= min_cells, , drop = FALSE]

logX <- log1p(count_f)
mu   <- Matrix::rowMeans(logX)
mu2  <- Matrix::rowMeans(logX^2)
v    <- pmax(mu2 - mu^2, 0)

n_hvg_real <- min(n_hvg, length(v))
hvg_genes  <- rownames(logX)[order(v, decreasing = TRUE)[seq_len(n_hvg_real)]]

n_cells_use <- min(n_cells_use_cap, ncol(count_f))

rm(obj, obj_sub, countMatrix, detected, logX, mu, mu2, v, md, keep_cells); gc()

log_msg(">>> 预处理完成: preHTC+HTC 细胞=", ncol(count_f),
        ", 检测>=", min_cells, "细胞的基因=", nrow(count_f), ", HVG=", length(hvg_genes))

run_gene <- function(nm, gene) {
  gdir <- file.path(out_root, paste0(nm, "_", gene, "过表达"))
  dir.create(gdir, showWarnings = FALSE, recursive = TRUE)
  res_file <- file.path(gdir, paste0(gene, "_preHTC_HTC_scTenifoldPerturb_OE_result.rds"))

  if (file.exists(res_file) && !force_rerun) {
    log_msg("SKIP ", nm, " ", gene, " (结果已存在)")
    return(invisible(FALSE))
  }

  t0 <- proc.time()
  log_msg("START ", nm, " ", gene)

  ok <- tryCatch({
    g_fixed <- fix_gene_symbol(gene, rownames(count_f))
    if (is.na(g_fixed)) stop("基因未通过 min_cells 过滤: ", gene)
    if (g_fixed != gene) log_msg("    基因名修正: ", gene, " -> ", g_fixed)

    count_small <- count_f[unique(c(g_fixed, hvg_genes)), , drop = FALSE]
    log_msg("    矩阵: ", nrow(count_small), " 基因 x ", ncol(count_small), " 细胞, 抽样 ", n_cells_use, " 细胞")

    res <- scTenifoldPerturb(
      countMatrix = count_small,
      gene        = g_fixed,
      mode        = "OE",
      perturb_factor = 2,
      qc          = FALSE,
      nc_lambda   = 0,
      nc_nNet     = 10,
      nc_nCells   = n_cells_use,
      nc_nComp    = 3,
      nc_scaleScores = TRUE,
      nc_symmetric   = FALSE,
      nc_q        = 0.9,
      td_K        = 3,
      td_maxIter  = 1000,
      td_maxError = 1e-5,
      td_nDecimal = 3,
      ma_nDim     = 2,
      seed        = 123,
      auto_adjust = FALSE,
      nCores      = nCores_use,
      safe_acceleration = TRUE,
      fast_pc     = FALSE,
      memory_efficient = TRUE
    )
    saveRDS(res, file = res_file)

    dr <- res$diffRegulation
    if (is.null(dr)) stop("Result object missing diffRegulation table.")

    if (!("gene" %in% colnames(dr))) {
      dr <- tibble::rownames_to_column(as.data.frame(dr), var = "gene")
    } else {
      dr <- as.data.frame(dr)
    }

    cn <- colnames(dr)
    p_col  <- cn[grep("padj|p\\.adj|FDR|adj\\.P\\.Val|adj_p|qval|q\\.val", cn, ignore.case = TRUE)][1]
    z_col  <- cn[grep("^Z$|^z$|zscore|z_score", cn, ignore.case = TRUE)][1]
    fc_col <- cn[grep("FC|foldChange|log2FC|logFC", cn, ignore.case = TRUE)][1]

    if (is.na(p_col) || is.na(z_col)) {
      stop("Cannot detect p-adj or Z column in diffRegulation. Columns: ", paste(cn, collapse = ", "))
    }

    df <- dr %>%
      dplyr::rename(p.adj = dplyr::all_of(p_col),
                    Z     = dplyr::all_of(z_col))
    if (!is.na(fc_col) && !(fc_col %in% c(p_col, z_col))) {
      df <- df %>% dplyr::rename(FC = dplyr::all_of(fc_col))
    }

    df$p.adj <- as.numeric(df$p.adj)
    df$Z     <- as.numeric(df$Z)
    df$log_padj <- -log10(df$p.adj + .Machine$double.eps)
    p_threshold <- 0.05
    df$Significance <- ifelse(df$p.adj < p_threshold, "Significant", "Not Significant")

    write.table(df, file.path(gdir, "diffRegulation_full.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)

    top_genes <- df %>% dplyr::arrange(desc(abs(Z))) %>% head(20)
    p_bar <- ggplot(top_genes, aes(x = Z, y = reorder(gene, Z), fill = Z > 0)) +
      geom_col() +
      scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"), guide = "none") +
      labs(title = paste0("Top 20 regulated genes (OE: ", g_fixed, "; preHTC+HTC)"),
           x = "Z-score", y = NULL) +
      theme_classic()
    ggsave(file.path(gdir, "1_Barplot_Top20_Z.pdf"), p_bar, width = 6, height = 6)

    label_genes <- df %>% dplyr::filter(p.adj < p_threshold) %>% dplyr::arrange(p.adj) %>% head(15) %>% pull(gene)
    label_genes <- unique(c(label_genes, "MMP13", g_fixed))

    p_vol <- ggplot(df, aes(x = Z, y = log_padj)) +
      geom_point(aes(color = Significance), alpha = 0.6, size = 1.5) +
      scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "#E64B35")) +
      geom_hline(yintercept = -log10(p_threshold), linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      ggrepel::geom_text_repel(data = subset(df, gene %in% label_genes),
                               aes(label = gene), max.overlaps = 50, size = 3) +
      labs(title = paste0("Volcano (OE: ", g_fixed, "; preHTC+HTC)"),
           x = "Z-score", y = "-log10(FDR)") +
      theme_classic() + theme(legend.position = "top")
    ggsave(file.path(gdir, "2_Volcano_Plot.pdf"), p_vol, width = 7, height = 6)

    pie_data <- df %>% dplyr::count(Significance) %>%
      dplyr::mutate(prop = n / sum(n) * 100, lab = paste0(round(prop, 1), "%"))
    p_pie <- ggplot(pie_data, aes(x = "", y = prop, fill = Significance)) +
      geom_bar(stat = "identity", width = 1) + coord_polar("y", start = 0) +
      geom_text(aes(label = lab), position = position_stack(vjust = 0.5)) +
      scale_fill_manual(values = c("Not Significant" = "#F0F0F0", "Significant" = "#E64B35")) +
      theme_void() +
      labs(title = "Significance proportion")
    ggsave(file.path(gdir, "3_Pie_Significant.pdf"), p_pie, width = 5, height = 5)

    if ("MMP13" %in% df$gene) {
      mmp_dat <- df %>% dplyr::filter(gene == "MMP13")
      p_mmp <- ggplot(mmp_dat, aes(x = Z, y = gene)) +
        geom_segment(aes(x = 0, xend = Z, y = gene, yend = gene),
                     linewidth = 1.2, color = "grey") +
        geom_point(size = 4, color = "#E64B35") +
        geom_text(aes(label = round(Z, 3)), vjust = -1.2) +
        labs(title = "MMP13 status under OE", x = "Z-score", y = NULL) +
        theme_minimal()
      ggsave(file.path(gdir, "4_MMP13_Lollipop.pdf"), p_mmp, width = 5, height = 3)
    }

    top_n_net <- 30
    net_data <- df %>%
      dplyr::filter(gene != g_fixed) %>%
      dplyr::arrange(desc(abs(Z))) %>%
      head(top_n_net)

    pdf(file.path(gdir, paste0("5_Network_OE_", g_fixed, "_preHTC_HTC.pdf")), width = 10, height = 8)
    tryCatch({
      edges <- data.frame(from = g_fixed, to = net_data$gene)
      g <- igraph::graph_from_data_frame(edges, directed = TRUE)

      V(g)$color <- "#E0E0E0"
      V(g)$size  <- 10
      V(g)$label.cex <- 0.8
      V(g)$frame.color <- "white"

      V(g)[g_fixed]$color <- "#3C5488"
      V(g)[g_fixed]$size  <- 22
      V(g)[g_fixed]$label.font <- 2
      V(g)[g_fixed]$label.color <- "white"

      node_names <- V(g)$name
      matches <- match(node_names, net_data$gene)
      z_values <- net_data$Z[matches]
      V(g)$color <- ifelse(is.na(z_values), "#3C5488", ifelse(z_values > 0, "#E64B35", "#4DBBD5"))

      l <- igraph::layout_as_star(g, center = g_fixed)
      plot(g,
           layout = l,
           edge.arrow.size = 0.5,
           edge.color = "grey80",
           vertex.label.color = ifelse(V(g)$name == g_fixed, "white", "black"),
           vertex.label.dist = 0,
           main = paste0("Virtual OE Network: ", g_fixed, " (preHTC+HTC)\nTop ", top_n_net, " regulated genes")
      )
      legend("bottomright",
             legend = c("OE Gene", "Up-regulated", "Down-regulated"),
             col = c("#3C5488", "#E64B35", "#4DBBD5"),
             pch = 19, pt.cex = 1.5, bty = "n")
    }, error = function(e) {
      log_msg("    网络图失败: ", e$message)
    })
    dev.off()

    log_msg("    表: ", nrow(df), " 行, 显著(FDR<0.05)=", sum(df$p.adj < 0.05),
            ", ", g_fixed, " 的 Z=", round(df$Z[df$gene == g_fixed], 3))
    TRUE
  }, error = function(e) {
    log_msg("ERROR ", nm, " ", gene, ": ", conditionMessage(e))
    FALSE
  })

  elapsed <- (proc.time() - t0)[3] / 60
  log_msg(if (isTRUE(ok)) "PASS" else "FAIL", " ", nm, " ", gene, " 耗时min=", round(elapsed, 1))
  while (!is.null(dev.list())) dev.off()
  gc()
  invisible(isTRUE(ok))
}

results <- sapply(names(gene_map), function(nm) run_gene(nm, unname(gene_map[nm])))

log_msg("=== 全部结束: 成功 ", sum(results), "/", length(results), " ===")
if (any(!results)) log_msg("失败: ", paste(names(gene_map)[!results], collapse = ", "))
