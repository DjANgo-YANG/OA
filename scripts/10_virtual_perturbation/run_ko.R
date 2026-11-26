rm(list = ls())
gc()

panel <- c(
  "AGL", "NASP", "LHX9", "COPB2", "PAM", "VCAN", "ELF3", "SDC2",
  "ME2", "PGK1", "HDLBP", "MDH1", "DCN", "PKM", "LDHA", "SLC35A3",
  "PLOD2", "CD44", "ENO1", "PPIA", "AK3", "TGFBI", "COG2"
)
args <- commandArgs(trailingOnly = TRUE)
plan_only <- "--plan" %in% args
selection <- args[args != "--plan"]
if (length(selection) != 1L || !selection %in% c("--all", panel)) {
  stop("Specify --all or one target from: ", paste(panel, collapse = ", "),
       ". Use --plan to list output paths without running the analysis.", call. = FALSE)
}

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
options(stringsAsFactors = FALSE)

target_genes <- if (selection == "--all") panel else selection
ko_root <- Sys.getenv("OAGLY_KO_ROOT", unset = project_path("results_root", "virtual_ko"))
target_dirs <- file.path(ko_root, sprintf("%02d_%s敲除", match(target_genes, panel), target_genes))
if (plan_only) {
  print(data.frame(target = target_genes, output_dir = target_dirs), row.names = FALSE)
  quit(save = "no")
}

in_rds <- project_path("results_root", "glycolysis_scoring/seurat_with_LeadingEdge_scores.rds")
if (!file.exists(in_rds)) stop("Input Seurat RDS was not found: ", in_rds)
existing <- target_dirs[file.exists(target_dirs)]
if (length(existing)) {
  stop("Output already exists. Choose a fresh OAGLY_KO_ROOT or run only missing targets: ",
       paste(existing, collapse = ", "))
}
dir.create(ko_root, recursive = TRUE, showWarnings = FALSE)

target_celltypes <- c("preHTC", "HTC")
min_cells <- 25
n_hvg <- 10000
n_cells_use_cap <- 800
nCores_use <- max(1, parallel::detectCores() - 1)

pkgs <- c("Matrix", "scTenifoldKnk", "ggplot2", "dplyr", "ggrepel", "igraph", "tibble", "stringr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing missing package: ", p)
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages(lapply(pkgs, library, character.only = TRUE))

if (!requireNamespace("SeuratObject", quietly = TRUE)) {
  install.packages("SeuratObject", repos = "https://cloud.r-project.org")
}

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
  if (inherits(x, "SingleCellExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) stop("Need SummarizedExperiment.")
    return(SummarizedExperiment::assay(x, "counts"))
  }
  stop("Unsupported object type for counts.")
}

fix_gene_symbol <- function(g, all_genes) {
  if (g %in% all_genes) return(g)
  tries <- unique(c(g, toupper(g), tolower(g), stringr::str_to_title(g)))
  hit <- tries[tries %in% all_genes]
  if (length(hit) > 0) return(hit[1])

  if (toupper(g) %in% all_genes) return(toupper(g))
  if (paste0(toupper(substr(g,1,1)), tolower(substr(g,2,nchar(g)))) %in% all_genes) {
    return(paste0(toupper(substr(g,1,1)), tolower(substr(g,2,nchar(g)))))
  }
  return(NA_character_)
}

message("Reading Seurat object: ", in_rds)
obj <- readRDS(in_rds)
if (!inherits(obj, "Seurat")) stop("Input RDS is not a Seurat object.")
ct_col <- pick_meta_col(obj, c("cell_type", "celltype", "CellType", "cellType", "annot", "annotation", "seurat_clusters"))
if (is.na(ct_col)) stop("Cannot find the cell type column in the Seurat metadata.")
md <- obj@meta.data
keep_cells <- rownames(md)[as.character(md[[ct_col]]) %in% target_celltypes]
if (length(keep_cells) < 50L) {
  stop("Too few cells after filtering preHTC/HTC (n=", length(keep_cells), "). Check cell type labels.")
}

counts <- get_counts(subset(obj, cells = keep_cells))
counts <- counts[Matrix::rowSums(counts != 0) >= min_cells, , drop = FALSE]
target_symbols <- vapply(target_genes, fix_gene_symbol, character(1), all_genes = rownames(counts))
if (anyNA(target_symbols)) {
  stop("Targets absent or detected in fewer than ", min_cells, " cells: ",
       paste(target_genes[is.na(target_symbols)], collapse = ", "))
}

log_counts <- log1p(counts)
variance <- pmax(Matrix::rowMeans(log_counts^2) - Matrix::rowMeans(log_counts)^2, 0)
hvg_genes <- rownames(counts)[head(order(variance, decreasing = TRUE), n_hvg)]
rm(obj, log_counts)
n_cells_use <- min(n_cells_use_cap, ncol(counts))

run_ko <- function(gKO, final_dir) {
  out_dir <- tempfile(pattern = ".ko-", tmpdir = ko_root)
  if (!dir.create(out_dir)) stop("Cannot create temporary KO output directory: ", out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  count_small <- counts[unique(c(gKO, hvg_genes)), , drop = FALSE]
  message("Running KO for ", gKO, " with ", nrow(count_small), " input genes.")
  set.seed(123)
  res <- scTenifoldKnk(
    countMatrix = count_small,
    gKO = gKO,
    qc = FALSE,
    nc_nNet = 10,
    nc_nCells = n_cells_use,
    nc_nComp = 3,
    nc_q = 0.9,
    td_K = 3,
    ma_nDim = 2,
    nCores = nCores_use
  )
  res_file <- file.path(out_dir, paste0(gKO, "_preHTC_HTC_scTenifoldKnk_result.rds"))
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

  write.table(df, file.path(out_dir, "diffRegulation_full.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  top_genes <- df %>% dplyr::arrange(desc(abs(Z))) %>% head(20)
  p_bar <- ggplot(top_genes, aes(x = Z, y = reorder(gene, Z), fill = Z > 0)) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"), guide = "none") +
    labs(title = paste0("Top 20 regulated genes (KO: ", gKO, "; preHTC+HTC)"),
         x = "Z-score", y = NULL) +
    theme_classic()
  ggsave(file.path(out_dir, "1_Barplot_Top20_Z.pdf"), p_bar, width = 6, height = 6)

  label_genes <- df %>% dplyr::filter(p.adj < p_threshold) %>% dplyr::arrange(p.adj) %>% head(15) %>% pull(gene)
  label_genes <- unique(c(label_genes, "MMP13", gKO))

  p_vol <- ggplot(df, aes(x = Z, y = log_padj)) +
    geom_point(aes(color = Significance), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "#E64B35")) +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    ggrepel::geom_text_repel(data = subset(df, gene %in% label_genes),
                             aes(label = gene), max.overlaps = 50, size = 3) +
    labs(title = paste0("Volcano (KO: ", gKO, "; preHTC+HTC)"),
         x = "Z-score", y = "-log10(FDR)") +
    theme_classic() + theme(legend.position = "top")
  ggsave(file.path(out_dir, "2_Volcano_Plot.pdf"), p_vol, width = 7, height = 6)

  pie_data <- df %>% dplyr::count(Significance) %>%
    dplyr::mutate(prop = n / sum(n) * 100, lab = paste0(round(prop, 1), "%"))
  p_pie <- ggplot(pie_data, aes(x = "", y = prop, fill = Significance)) +
    geom_bar(stat = "identity", width = 1) + coord_polar("y", start = 0) +
    geom_text(aes(label = lab), position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = c("Not Significant" = "#F0F0F0", "Significant" = "#E64B35")) +
    theme_void() +
    labs(title = "Significance proportion")
  ggsave(file.path(out_dir, "3_Pie_Significant.pdf"), p_pie, width = 5, height = 5)

  if ("MMP13" %in% df$gene) {
    mmp_dat <- df %>% dplyr::filter(gene == "MMP13")
    p_mmp <- ggplot(mmp_dat, aes(x = Z, y = gene)) +
      geom_segment(aes(x = 0, xend = Z, y = gene, yend = gene),
                   linewidth = 1.2, color = "grey") +
      geom_point(size = 4, color = "#E64B35") +
      geom_text(aes(label = round(Z, 3)), vjust = -1.2) +
      labs(title = "MMP13 status under KO", x = "Z-score", y = NULL) +
      theme_minimal()
    ggsave(file.path(out_dir, "4_MMP13_Lollipop.pdf"), p_mmp, width = 5, height = 3)
  }

  message(">>> Drawing KO-centered network (local igraph) ...")
  top_n_net <- 30
  net_data <- df %>%
    dplyr::filter(gene != gKO) %>%
    dplyr::arrange(desc(abs(Z))) %>%
    head(top_n_net)

  pdf(file.path(out_dir, paste0("5_Network_KO_", gKO, "_preHTC_HTC.pdf")), width = 10, height = 8)
  tryCatch({
    edges <- data.frame(from = gKO, to = net_data$gene)
    g <- igraph::graph_from_data_frame(edges, directed = TRUE)

    V(g)$color <- "#E0E0E0"
    V(g)$size  <- 10
    V(g)$label.cex <- 0.8
    V(g)$frame.color <- "white"

    V(g)[gKO]$color <- "#3C5488"
    V(g)[gKO]$size  <- 22
    V(g)[gKO]$label.font <- 2
    V(g)[gKO]$label.color <- "white"

    node_names <- V(g)$name
    matches <- match(node_names, net_data$gene)
    z_values <- net_data$Z[matches]
    V(g)$color <- ifelse(is.na(z_values), "#3C5488", ifelse(z_values > 0, "#E64B35", "#4DBBD5"))

    l <- igraph::layout_as_star(g, center = gKO)
    plot(g,
         layout = l,
         edge.arrow.size = 0.5,
         edge.color = "grey80",
         vertex.label.color = ifelse(V(g)$name == gKO, "white", "black"),
         vertex.label.dist = 0,
         main = paste0("Virtual KO Network: ", gKO, " (preHTC+HTC)\nTop ", top_n_net, " regulated genes")
    )
    legend("bottomright",
           legend = c("KO Gene", "Up-regulated", "Down-regulated"),
           col = c("#3C5488", "#E64B35", "#4DBBD5"),
           pch = 19, pt.cex = 1.5, bty = "n")
  }, error = function(e) {
    message("Network plot failed: ", e$message)
  })
  dev.off()

  if (file.exists(final_dir) || !file.rename(out_dir, final_dir)) {
    stop("Could not publish KO output: ", final_dir)
  }
  message("Saved KO results: ", final_dir)
}

for (i in seq_along(target_genes)) {
  run_ko(target_symbols[[i]], target_dirs[[i]])
}
