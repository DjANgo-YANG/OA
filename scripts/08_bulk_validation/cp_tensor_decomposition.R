source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))

options(stringsAsFactors = FALSE)

WORKDIR    <- project_path("data_root", "bulk")
COUNT_FILE <- file.path(WORKDIR, "GSE254844_processed_counts.txt")
OUT_ROOT   <- file.path(WORKDIR, "Tensor_CP_noSampleAxis")

dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)
setwd(OUT_ROOT)

N_TOP      <- 3000
GENE_PICK  <- "var"
RANK       <- 4
TOPK       <- 200

DO_RANK_SELECTION <- FALSE
R_try <- 1:7
SEEDS <- c(1)

MIN_SIZE   <- 10
MAX_SIZE   <- 500
TOP_SHOW   <- 30

TOP_PLOT     <- 6
PLOT_FDR_CUT <- 0.25

COL_DARK  <- "#08519c"
COL_LIGHT <- "#08519c"

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(dplyr)
  library(tibble)
  library(rTensor)
  library(matrixStats)
  library(msigdbr)
  library(fgsea)
  library(enrichplot)
  library(ggplot2)
  library(tidyr)
})

OA_COUNTS <- read.table(
  COUNT_FILE,
  header       = TRUE,
  sep          = "",
  quote        = "",
  comment.char = "",
  check.names  = FALSE,
  row.names    = 1
)
write.csv(OA_COUNTS, "QC_OA_COUNTS_raw.csv")
sample <- colnames(OA_COUNTS)

condition <- ifelse(grepl("^Control\\.", sample), "HC",
                    ifelse(grepl("^(WB|NWB)\\.", sample), "OA", NA))

layer <- sub("^(Control|NWB|WB)\\.([A-Z]{2}).*$", "\\2", sample)
layer[!layer %in% c("AS","SZ","MZ","DZ")] <- NA

if (any(is.na(condition)) | any(is.na(layer))) {
  bad <- sample[is.na(condition) | is.na(layer)]
  stop("解析失败样本示例：", paste(head(bad, 10), collapse = ", "))
}

sample_info <- data.frame(sample_id = sample, condition = condition, layer = layer)
write.csv(sample_info, "sample_info_condition_layer.csv", row.names = FALSE)

cat("Sample counts by condition×layer:\n")
print(table(sample_info$condition, sample_info$layer))

cond_levels  <- c("HC","OA")
layer_levels <- c("AS","SZ","MZ","DZ")

for (cnd in cond_levels) for (lay in layer_levels) {
  if (sum(sample_info$condition == cnd & sample_info$layer == lay) == 0) {
    stop("缺少抽屉组合：", cnd, " × ", lay)
  }
}

expressed <- rowMeans(OA_COUNTS > 0) >= 0.2
OA_COUNTS_f <- OA_COUNTS[expressed, , drop = FALSE]
cat("Genes kept after filtering:", nrow(OA_COUNTS_f), "\n")

dge <- DGEList(counts = OA_COUNTS_f)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE, prior.count = 1)
stopifnot(all(colnames(logCPM) == sample_info$sample_id))

if (GENE_PICK == "var") {
  gene_score <- matrixStats::rowVars(logCPM)
  top_genes <- names(sort(gene_score, decreasing = TRUE))[1:min(N_TOP, length(gene_score))]
} else if (GENE_PICK == "t") {
  design <- model.matrix(~ 0 + factor(sample_info$condition, levels = c("HC","OA")))
  colnames(design) <- c("HC","OA")
  fit <- lmFit(logCPM, design)
  fit <- contrasts.fit(fit, makeContrasts(OA_vs_HC = OA - HC, levels = design))
  fit <- eBayes(fit, trend = TRUE)
  tt <- topTable(fit, coef = "OA_vs_HC", number = Inf, sort.by = "none")
  gene_score <- abs(tt$t)
  names(gene_score) <- rownames(tt)
  gene_score <- gene_score[rownames(logCPM)]
  top_genes <- names(sort(gene_score, decreasing = TRUE))[1:min(N_TOP, length(gene_score))]
} else {
  stop("GENE_PICK must be 'var' or 't'.")
}
writeLines(top_genes, "TopGenes_used_for_tensor.txt")

mean_arr_gene_cond_layer <- array(
  NA_real_,
  dim = c(length(top_genes), length(cond_levels), length(layer_levels)),
  dimnames = list(top_genes, cond_levels, layer_levels)
)

for (cnd in cond_levels) for (lay in layer_levels) {
  idx <- which(sample_info$condition == cnd & sample_info$layer == lay)
  mean_arr_gene_cond_layer[, cnd, lay] <- rowMeans(logCPM[top_genes, idx, drop = FALSE])
}

Y <- aperm(mean_arr_gene_cond_layer, c(2, 3, 1))
cat("Tensor Y dim:", paste(dim(Y), collapse = " x "), "\n")

Y_mat <- matrix(Y, nrow = length(cond_levels) * length(layer_levels), byrow = FALSE)
colnames(Y_mat) <- top_genes
Y_mat_centered <- scale(Y_mat, center = TRUE, scale = FALSE)

Y_centered <- array(
  Y_mat_centered,
  dim = c(length(cond_levels), length(layer_levels), length(top_genes)),
  dimnames = list(cond_levels, layer_levels, top_genes)
)

Y_tensor <- as.tensor(Y_centered)

frob_norm <- function(x) sqrt(sum(x * x))

reconstruct_cp_mat <- function(cp_fit, I, J, K) {
  U <- cp_fit$U
  lambda <- cp_fit$lambda
  R <- length(lambda)

  A <- U[[1]]
  B <- U[[2]]
  C <- U[[3]]

  Mhat <- matrix(0, nrow = I * J, ncol = K)
  for (r in 1:R) {
    u <- as.vector(outer(A[, r], B[, r]))
    Mhat <- Mhat + lambda[r] * (u %*% t(C[, r]))
  }
  Mhat
}

rel_error_cp <- function(Y_true_arr, cp_fit) {
  I <- dim(Y_true_arr)[1]
  J <- dim(Y_true_arr)[2]
  K <- dim(Y_true_arr)[3]
  M_true <- matrix(Y_true_arr, nrow = I * J, ncol = K)
  M_hat  <- reconstruct_cp_mat(cp_fit, I, J, K)
  frob_norm(M_true - M_hat) / frob_norm(M_true)
}

msig_h <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)

pathways_h <- split(msig_h$gene_symbol, msig_h$gs_name)
pathways_h <- lapply(pathways_h, function(gs) intersect(gs, top_genes))
pathways_h <- pathways_h[lengths(pathways_h) >= MIN_SIZE]

top_hallmark_from_loadings <- function(loadings_vec) {
  ranks <- loadings_vec
  names(ranks) <- names(loadings_vec)
  ranks <- ranks[!is.na(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)
  fg <- fgsea(pathways = pathways_h, stats = ranks, minSize = MIN_SIZE, maxSize = MAX_SIZE)
  fg <- fg[order(-abs(fg$NES), fg$padj), ]
  if (nrow(fg) == 0) return(NA_character_)
  fg$pathway[1]
}

write_fgsea <- function(fg, out_csv) {
  fg_out <- as.data.frame(fg)
  list_cols <- names(fg_out)[vapply(fg_out, is.list, logical(1))]
  if (length(list_cols) > 0) {
    for (cc in list_cols) {
      fg_out[[cc]] <- vapply(
        fg_out[[cc]],
        function(x) paste(x, collapse = ";"),
        FUN.VALUE = character(1)
      )
    }
  }
  write.csv(fg_out, out_csv, row.names = FALSE)
  fg_out
}

plot_factor <- function(mat, title, out_png) {
  df <- as.data.frame(mat) %>% rownames_to_column("level")

  df2 <- tidyr::pivot_longer(
    df,
    cols = -level,
    names_to = "comp",
    values_to = "weight"
  ) %>%
    dplyr::group_by(comp) %>%
    dplyr::mutate(fill_value = rank(weight, ties.method = "first")) %>%
    dplyr::ungroup()

  p <- ggplot(df2, aes(x = level, y = weight, fill = fill_value)) +
    geom_col(width = 0.8) +
    facet_wrap(~ comp, scales = "free_y") +
    scale_fill_gradient(low = COL_LIGHT, high = COL_DARK) +
    theme_bw() +
    labs(title = title, x = NULL, y = "weight") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )

  ggsave(out_png, p, width = 8, height = 3.8, dpi = 300)
}

plot_fgsea_bar <- function(fg_out, title, out_png, n = TOP_SHOW) {
  if (is.null(fg_out) || nrow(fg_out) == 0) return(invisible(NULL))

  df2 <- fg_out %>%
    dplyr::filter(!is.na(padj), !is.na(NES)) %>%
    dplyr::arrange(padj, desc(abs(NES))) %>%
    dplyr::slice_head(n = n)

  if (nrow(df2) == 0) return(invisible(NULL))

  df2 <- df2 %>%
    dplyr::mutate(
      pathway = factor(pathway, levels = rev(pathway)),
      fill_value = rank(abs(NES), ties.method = "first")
    )

  p <- ggplot(df2, aes(x = pathway, y = NES, fill = fill_value)) +
    geom_col(width = 0.8) +
    coord_flip() +
    scale_fill_gradient(low = COL_LIGHT, high = COL_DARK) +
    theme_bw() +
    labs(title = title, x = NULL, y = "NES") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )

  ggsave(out_png, p, width = 10, height = 7, dpi = 300)
}

if (DO_RANK_SELECTION) {

  rank_rows <- list()

  for (R0 in R_try) {
    for (sd in SEEDS) {
      set.seed(sd)
      cp_fit0 <- cp(Y_tensor, num_components = R0, max_iter = 500, tol = 1e-7)
      err0 <- rel_error_cp(Y_centered, cp_fit0)

      C0 <- cp_fit0$U[[3]]
      rownames(C0) <- top_genes
      colnames(C0) <- paste0("Comp", 1:R0)

      tops0 <- sapply(1:R0, function(k) {
        v <- C0[, k]
        names(v) <- rownames(C0)
        top_hallmark_from_loadings(v)
      })

      df0 <- data.frame(
        seed = sd,
        RANK = R0,
        RelError = err0,
        Component = paste0("Comp", 1:R0),
        TopHallmark = as.character(tops0),
        stringsAsFactors = FALSE
      )
      rank_rows[[paste(R0, sd, sep = "_")]] <- df0

      cat("RANK=", R0, " seed=", sd, " RelError=", err0, "\n")
    }
  }

  rank_report <- bind_rows(rank_rows)
  write.csv(rank_report, "TensorRank_selection_report.csv", row.names = FALSE)

  curve <- rank_report %>%
    group_by(RANK) %>%
    summarise(RelError_min = min(RelError), RelError_median = median(RelError), .groups = "drop") %>%
    arrange(RANK)
  write.csv(curve, "rank_curve.csv", row.names = FALSE)

  suggest_elbow <- function(df, x_col = "RANK", y_col = "RelError_median") {
    x <- df[[x_col]]
    y <- df[[y_col]]

    x1 <- x[1]; y1 <- y[1]
    x2 <- x[length(x)]; y2 <- y[length(y)]

    num <- abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1)
    den <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
    d <- num / den

    x[which.max(d)]
  }

  elbow_R <- suggest_elbow(curve, "RANK", "RelError_median")
  writeLines(paste0("Elbow_suggested_RANK=", elbow_R), "elbow_suggestion.txt")
  cat("AUTO elbow suggested RANK =", elbow_R, "\n")

  p <- ggplot(curve, aes(x = RANK, y = RelError_median)) +
    geom_line() +
    geom_point() +
    theme_bw() +
    labs(title = "RANK selection curve (median RelError)", y = "RelError", x = "RANK") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave("rank_curve.png", p, width = 6, height = 4, dpi = 300)

  cat("\n[STOP] Rank selection done.\n",
      "Check: rank_curve.png / rank_curve.csv and TensorRank_selection_report.csv\n",
      "Then set DO_RANK_SELECTION <- FALSE and set RANK <- <your choice>, re-run.\n")

  stop("Stop after RANK selection. Set RANK manually and rerun.")
}

cp_fit <- cp(Y_tensor, num_components = RANK, max_iter = 800, tol = 1e-8)

A <- cp_fit$U[[1]]
B <- cp_fit$U[[2]]
C <- cp_fit$U[[3]]

rownames(A) <- cond_levels
rownames(B) <- layer_levels
rownames(C) <- top_genes
colnames(A) <- paste0("Comp", 1:RANK)
colnames(B) <- paste0("Comp", 1:RANK)
colnames(C) <- paste0("Comp", 1:RANK)

write.csv(A, "Tensor_CP_condition_factors.csv")
write.csv(B, "Tensor_CP_layer_factors.csv")
write.csv(C, "Tensor_CP_gene_loadings.csv")

plot_factor(A, "CP factors: condition", "Tensor_CP_condition_factors.png")
plot_factor(B, "CP factors: layer", "Tensor_CP_layer_factors.png")

dir.create("TopGenes", showWarnings = FALSE)
for (r in 1:RANK) {
  w <- C[, r]
  genes_r <- names(sort(abs(w), decreasing = TRUE))[1:TOPK]
  writeLines(genes_r, file.path("TopGenes", paste0("Tensor_Comp", r, "_Top", TOPK, ".txt")))
}

dir.create("Hallmark_GSEA", showWarnings = FALSE)

for (r in 1:RANK) {
  comp_name <- paste0("Comp", r)
  ranks <- C[, r]
  names(ranks) <- rownames(C)
  ranks <- sort(ranks, decreasing = TRUE)

  fg <- fgsea(pathways = pathways_h, stats = ranks, minSize = MIN_SIZE, maxSize = MAX_SIZE)
  fg <- fg[order(fg$padj, -abs(fg$NES)), ]

  fg_out <- write_fgsea(fg, file.path("Hallmark_GSEA", paste0(comp_name, "_Hallmark_GSEA.csv")))

  top_tbl <- fg_out %>%
    dplyr::filter(!is.na(padj), !is.na(NES)) %>%
    dplyr::arrange(padj, desc(abs(NES))) %>%
    dplyr::slice_head(n = TOP_SHOW)

  write.csv(
    top_tbl,
    file.path("Hallmark_GSEA", paste0(comp_name, "_Hallmark_TOP", TOP_SHOW, ".csv")),
    row.names = FALSE
  )

  plot_fgsea_bar(
    fg_out,
    title = paste0(comp_name, " Hallmark GSEA"),
    out_png = file.path("Hallmark_GSEA", paste0(comp_name, "_Hallmark_TOP", TOP_SHOW, "_barplot.png")),
    n = TOP_SHOW
  )

  pick <- fg_out %>%
    dplyr::filter(!is.na(padj), !is.na(NES), padj < PLOT_FDR_CUT) %>%
    dplyr::arrange(padj, desc(abs(NES))) %>%
    dplyr::slice_head(n = TOP_PLOT) %>%
    dplyr::pull(pathway)

  pick <- intersect(pick, names(pathways_h))

  if (length(pick) > 0) {
    pdf(
      file.path("Hallmark_GSEA", paste0(comp_name, "_AutoTop", TOP_PLOT, "_EnrichmentPlots.pdf")),
      width = 7,
      height = 5
    )
    for (pw in pick) {
      print(
        plotEnrichment(pathways_h[[pw]], ranks) +
          labs(title = paste0(comp_name, " | ", pw))
      )
    }
    dev.off()
  }
}

cat("\nDONE.\nOutputs in: ", OUT_ROOT, "\n",
    "Key files:\n",
    " - rank_curve.png / rank_curve.csv (if rank selection run)\n",
    " - Tensor_CP_condition_factors.csv (+png)\n",
    " - Tensor_CP_layer_factors.csv (+png)\n",
    " - Hallmark_GSEA/\n")
