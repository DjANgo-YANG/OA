rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
options(stringsAsFactors = FALSE)

INPUT_RDS <- project_path("data_root", "auxiliary/reference/glycolysis_scoring/seurat_with_LeadingEdge_scores.rds")
OUT_DIR   <- project_path("results_root", "exploratory_perturbation")

CELLTYPES <- c("HTC", "preHTC")
CONDITION <- "OA"

GLYCOLYSIS_GENES <- c(
  "HK1", "HK2", "HK3",
  "GPI", "PFKL", "PFKM", "PFKP",
  "ALDOA", "ALDOB", "ALDOC",
  "TPI1",
  "GAPDH",
  "PGK1",
  "PGAM1", "PGAM2", "PGAM4",
  "ENO1", "ENO2", "ENO3",
  "PKM", "PKLR",
  "LDHA", "LDHB", "LDHC",
  "PDK1", "PDK2", "PDK3", "PDK4",
  "PDHA1", "PDHA2", "PDHB",
  "SLC2A1", "SLC2A3", "SLC2A4",
  "PFKFB1", "PFKFB2", "PFKFB3", "PFKFB4",
  "G6PC", "G6PC2", "G6PC3"
)

set.seed(1234)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

message("Reading Seurat object ...")
seu <- readRDS(INPUT_RDS)

if (!"condition" %in% colnames(seu@meta.data)) {
  seu@meta.data$condition <- ifelse(grepl("^HC", seu@meta.data$orig.ident), "HC", "OA")
  seu@meta.data$condition <- factor(seu@meta.data$condition, levels = c("HC", "OA"))
}

message("\n=== Starting Analysis ===\n")

all_results <- list()
summary_list <- list()

for (ct in CELLTYPES) {
  cat("\n========================================\n")
  cat("Cell Type:", ct, "\n")
  cat("========================================\n")

  ct_dir <- file.path(OUT_DIR, ct)
  if (!dir.exists(ct_dir)) dir.create(ct_dir, recursive = TRUE)

  cells_ct <- rownames(seu@meta.data)[seu@meta.data$cell_type == ct]
  cat("Cells in", ct, ":", length(cells_ct), "\n")

  cells_oa <- cells_ct[seu@meta.data[cells_ct, "condition"] == "OA"]
  cells_hc <- cells_ct[seu@meta.data[cells_ct, "condition"] == "HC"]

  cat("OA cells:", length(cells_oa), "\n")
  cat("HC cells:", length(cells_hc), "\n")

  expr <- seu@assays$RNA@data

  expr_oa <- expr[, cells_oa, drop = FALSE]
  expr_hc <- expr[, cells_hc, drop = FALSE]

  glycolysis_present <- intersect(GLYCOLYSIS_GENES, rownames(expr_oa))
  cat("Glycolysis genes present:", length(glycolysis_present), "\n")

  glycolysis_stats <- data.frame(
    gene = glycolysis_present,
    mean_expr_oa = rowMeans(expr_oa[glycolysis_present, , drop = FALSE]),
    mean_expr_hc = rowMeans(expr_hc[glycolysis_present, , drop = FALSE]),
    pct_expr_oa = rowMeans(expr_oa[glycolysis_present, , drop = FALSE] > 0),
    pct_expr_hc = rowMeans(expr_hc[glycolysis_present, , drop = FALSE] > 0),
    logFC = rowMeans(expr_oa[glycolysis_present, , drop = FALSE]) - rowMeans(expr_hc[glycolysis_present, , drop = FALSE])
  )

  glycolysis_stats <- glycolysis_stats[order(glycolysis_stats$mean_expr_oa, decreasing = TRUE), ]

  write.csv(glycolysis_stats, file.path(ct_dir, "glycolysis_genes_OA_vs_HC.csv"), row.names = FALSE)

  core_threshold <- quantile(glycolysis_stats$mean_expr_oa, 0.75)
  core_genes <- glycolysis_stats[glycolysis_stats$mean_expr_oa >= core_threshold, "gene"]
  cat("Core glycolysis genes (high expr in OA):", length(core_genes), "\n")
  writeLines(as.character(core_genes), file.path(ct_dir, "core_glycolysis_genes.txt"))

  glycolysis_genes_in_data <- intersect(GLYCOLYSIS_GENES, rownames(expr_oa))
  glycolysis_score_baseline <- colMeans(expr_oa[glycolysis_genes_in_data, , drop = FALSE], na.rm = TRUE)

  cat("\n--- Virtual Knockout Analysis ---\n")
  ko_results <- data.frame(
    gene = character(),
    cell_type = character(),
    operation = character(),
    glycolysis_score_orig = numeric(),
    glycolysis_score_after = numeric(),
    glycolysis_delta = numeric(),
    stringsAsFactors = FALSE
  )

  for (gene in core_genes) {
    if (!(gene %in% rownames(expr_oa))) next

    expr_ko <- expr_oa
    expr_ko[gene, ] <- 0

    glycolysis_score_ko <- colMeans(expr_ko[glycolysis_genes_in_data, , drop = FALSE], na.rm = TRUE)

    ko_row <- data.frame(
      gene = gene,
      cell_type = ct,
      operation = "KO",
      glycolysis_score_orig = mean(glycolysis_score_baseline),
      glycolysis_score_after = mean(glycolysis_score_ko),
      glycolysis_delta = mean(glycolysis_score_ko) - mean(glycolysis_score_baseline),
      stringsAsFactors = FALSE
    )
    ko_results <- rbind(ko_results, ko_row)
  }

  write.csv(ko_results, file.path(ct_dir, "virtual_KO_results.csv"), row.names = FALSE)

  cat("\n--- Virtual Overexpression Analysis ---\n")
  oe_results <- data.frame(
    gene = character(),
    cell_type = character(),
    operation = character(),
    glycolysis_score_orig = numeric(),
    glycolysis_score_after = numeric(),
    glycolysis_delta = numeric(),
    stringsAsFactors = FALSE
  )

  for (gene in core_genes) {
    if (!(gene %in% rownames(expr_oa))) next

    expr_oe <- expr_oa
    expr_oe[gene, ] <- expr_oe[gene, ] * 2

    glycolysis_score_oe <- colMeans(expr_oe[glycolysis_genes_in_data, , drop = FALSE], na.rm = TRUE)

    oe_row <- data.frame(
      gene = gene,
      cell_type = ct,
      operation = "OE",
      glycolysis_score_orig = mean(glycolysis_score_baseline),
      glycolysis_score_after = mean(glycolysis_score_oe),
      glycolysis_delta = mean(glycolysis_score_oe) - mean(glycolysis_score_baseline),
      stringsAsFactors = FALSE
    )
    oe_results <- rbind(oe_results, oe_row)
  }

  write.csv(oe_results, file.path(ct_dir, "virtual_OE_results.csv"), row.names = FALSE)

  summary_stats <- data.frame(
    cell_type = ct,
    n_cells_oa = length(cells_oa),
    n_cells_hc = length(cells_hc),
    n_glycolysis_genes = length(glycolysis_present),
    n_core_genes = length(core_genes),
    mean_glycolysis_score_oa = mean(glycolysis_score_baseline),
    ko_mean_delta = mean(ko_results$glycolysis_delta),
    oe_mean_delta = mean(oe_results$glycolysis_delta),
    stringsAsFactors = FALSE
  )

  summary_list[[ct]] <- summary_stats
  write.csv(summary_stats, file.path(ct_dir, "summary_statistics.csv"), row.names = FALSE)

  cat("\n", ct, " Analysis Complete\n")
  cat("Core genes:", length(core_genes), "\n")
  cat("KO mean delta:", round(mean(ko_results$glycolysis_delta), 4), "\n")
  cat("OE mean delta:", round(mean(oe_results$glycolysis_delta), 4), "\n")
}

summary_df <- do.call(rbind, summary_list)
write.csv(summary_df, file.path(OUT_DIR, "overall_summary.csv"), row.names = FALSE)

message("\n========================================")
message("All analyses complete!")
message("Output directory: ", OUT_DIR)
message("========================================\n")

message("Summary of findings:\n")
print(summary_df)
