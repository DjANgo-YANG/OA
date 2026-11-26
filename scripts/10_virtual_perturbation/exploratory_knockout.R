source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))

in_rds  <- project_path("data_root", "auxiliary/oa/scored_subset.rds")
out_dir <- project_path("results_root", "exploratory_knockout")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

gKO <- "HIF1A"

pkgs <- c("Matrix", "scTenifoldKnk")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(Matrix)
library(scTenifoldKnk)

if (!requireNamespace("SeuratObject", quietly = TRUE)) {

}
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {

}
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {

}

obj <- readRDS(in_rds)

get_counts <- function(x) {
  if (inherits(x, "dgCMatrix") || is.matrix(x)) return(x)

  if (inherits(x, "Seurat")) {
    m <- SeuratObject::GetAssayData(x, assay = "RNA", slot = "counts")
    return(m)
  }

  if (inherits(x, "SingleCellExperiment")) {
    m <- SummarizedExperiment::assay(x, "counts")
    return(m)
  }

  stop("不认识的RDS对象类型：既不是矩阵，也不是Seurat/SCE。请提供 str(obj) 前几行。")
}

countMatrix <- get_counts(obj)

if (is.null(rownames(countMatrix)) || is.null(colnames(countMatrix))) {
  stop("countMatrix 必须有 rownames(基因名) 和 colnames(细胞名)。")
}

if (!gKO %in% rownames(countMatrix)) {
  message("找不到基因：", gKO)
  gKO2 <- ifelse(gKO == toupper(gKO),
                 paste0(tolower(substr(gKO, 1, 1)), tolower(substr(gKO, 2, nchar(gKO)))),
                 toupper(gKO))
  message("尝试变体：", gKO2)
  if (gKO2 %in% rownames(countMatrix)) {
    gKO <- gKO2
    message("已改用：", gKO)
  } else {
    stop("ETV1/Etv1 不在表达矩阵行名里。请确认物种基因符号或先做gene symbol转换。")
  }
}

set.seed(1)

min_cells <- 25
n_hvg     <- 5000
keep_always <- unique(c(
  gKO
))
keep_always <- keep_always[keep_always %in% rownames(countMatrix)]

detected <- Matrix::rowSums(countMatrix != 0)
count_f  <- countMatrix[detected >= min_cells, , drop = FALSE]

logX <- log1p(count_f)
mu  <- Matrix::rowMeans(logX)
mu2 <- Matrix::rowMeans(logX^2)
v   <- pmax(mu2 - mu^2, 0)

n_hvg <- min(n_hvg, length(v))
hvg_genes <- rownames(logX)[order(v, decreasing = TRUE)[seq_len(n_hvg)]]

keep_genes <- unique(c(keep_always, hvg_genes))
count_small <- count_f[keep_genes, , drop = FALSE]

message("After gene selection: genes=", nrow(count_small), ", cells=", ncol(count_small))
if (nrow(count_small) > 8000) {
  message("警告：基因仍>8000，建议把 n_hvg 调小（如 3000–5000）以避免内存问题。")
}

n_cells <- ncol(count_small)
nc_nCells_use <- min(500, n_cells)

res <- scTenifoldKnk(
  countMatrix = count_small,
  gKO         = gKO,
  qc          = FALSE,
  nc_nNet     = 10,
  nc_nCells   = nc_nCells_use,
  nc_nComp    = 3,
  nc_q        = 0.9,
  td_K        = 3,
  ma_nDim     = 2,
  nCores      = max(1, parallel::detectCores() - 1)
)

saveRDS(res, file = sprintf("ETV1_scTenifoldKnk_result_%s.rds", gKO))

dr <- res$diffRegulation
write.csv(dr, file = sprintf("ETV1_virtualKO_diffRegulation_%s.csv", gKO), row.names = FALSE)

writeLines(paste("diffRegulation columns:", paste(colnames(dr), collapse = ", ")),
           con = "diffRegulation_columns.txt")

ord <- NULL
if ("padj" %in% colnames(dr)) ord <- order(dr$padj, na.last = TRUE)
if (is.null(ord) && "FDR" %in% colnames(dr)) ord <- order(dr$FDR, na.last = TRUE)
if (is.null(ord) && "pvalue" %in% colnames(dr)) ord <- order(dr$pvalue, na.last = TRUE)
if (is.null(ord)) {
  num_cols <- which(sapply(dr, is.numeric))
  if (length(num_cols) > 0) ord <- order(dr[[num_cols[1]]], decreasing = TRUE, na.last = TRUE)
}
if (!is.null(ord)) {
  top <- dr[ord, , drop = FALSE]
  write.csv(head(top, 200), file = sprintf("ETV1_virtualKO_top200_%s.csv", gKO), row.names = FALSE)
}

pdf(sprintf("ETV1_virtualKO_network_%s.pdf", gKO), width = 10, height = 8)
try({
  plotKO(res, gKO = gKO, q = 0.99, annotate = FALSE, fdrThreshold = 0.05)
}, silent = TRUE)
dev.off()

message("DONE. Outputs written to: ", out_dir)
