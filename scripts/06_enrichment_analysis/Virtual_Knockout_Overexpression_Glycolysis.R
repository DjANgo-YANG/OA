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

OA_PATHWAYS <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_HYPOXIA",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_APOPTOSIS",
  "GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY",
  "GOBP_COLLAGEN_FIBRIL_ORGANIZATION",
  "GOBP_CHONDROCYTE_DIFFERENTIATION"
)

set.seed(1234)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

pkgs_cran <- c("ggplot2", "dplyr", "tibble", "patchwork")
pkgs_bioc <- c("Seurat", "msigdbr", "clusterProfiler", "enrichplot")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in pkgs_cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}
for (p in pkgs_bioc) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, update = FALSE, ask = FALSE)
  library(p, character.only = TRUE)
}

message("Reading Seurat object ...")
seu <- readRDS(INPUT_RDS)

if (!"condition" %in% colnames(seu@meta.data)) {
  seu@meta.data$condition <- ifelse(grepl("^HC", seu@meta.data$orig.ident), "HC", "OA")
  seu@meta.data$condition <- factor(seu@meta.data$condition, levels = c("HC", "OA"))
}

message("Loading gene sets ...")
m_H  <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
m_C5 <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
m_t2g <- rbind(m_H, m_C5)

hallmark_glycolysis <- m_H[m_H$gs_name == "HALLMARK_GLYCOLYSIS", "gene_symbol"]
ALL_GLYCOLYSIS_GENES <- unique(c(GLYCOLYSIS_GENES, hallmark_glycolysis))
message("Total glycolysis genes: ", length(ALL_GLYCOLYSIS_GENES))

oa_pathway_genes <- list()
for (pathway in OA_PATHWAYS) {
  genes <- m_t2g[m_t2g$gs_name == pathway, "gene_symbol"]
  oa_pathway_genes[[pathway]] <- genes
  message(pathway, ": ", length(genes), " genes")
}

message("\n=== Identifying core glycolysis genes ===\n")

for (ct in CELLTYPES) {
  ct_dir <- file.path(OUT_DIR, ct)
  if (!dir.exists(ct_dir)) dir.create(ct_dir, recursive = TRUE)

  cells_ct <- rownames(seu@meta.data)[seu@meta.data$cell_type == ct]
  seu_ct <- subset(seu, cells = cells_ct)

  cells_oa <- rownames(seu_ct@meta.data)[seu_ct@meta.data$condition == "OA"]
  seu_oa <- subset(seu_ct, cells = cells_oa)

  cells_hc <- rownames(seu_ct@meta.data)[seu_ct@meta.data$condition == "HC"]
  seu_hc <- subset(seu_ct, cells = cells_hc)

  DefaultAssay(seu_oa) <- "RNA"
  DefaultAssay(seu_hc) <- "RNA"

  expr_oa <- GetAssayData(seu_oa, slot = "data")
  expr_hc <- GetAssayData(seu_hc, slot = "data")

  glycolysis_present <- intersect(ALL_GLYCOLYSIS_GENES, rownames(expr_oa))
  message(ct, ": ", length(glycolysis_present), " glycolysis genes present")

  glycolysis_stats <- data.frame(
    gene = glycolysis_present,
    mean_expr_oa = rowMeans(expr_oa[glycolysis_present, ]),
    mean_expr_hc = rowMeans(expr_hc[glycolysis_present, ]),
    pct_expr_oa = rowMeans(expr_oa[glycolysis_present, ] > 0),
    pct_expr_hc = rowMeans(expr_hc[glycolysis_present, ] > 0),
    logFC = rowMeans(expr_oa[glycolysis_present, ]) - rowMeans(expr_hc[glycolysis_present, ])
  )

  glycolysis_stats <- glycolysis_stats[order(glycolysis_stats$mean_expr_oa, decreasing = TRUE), ]

  write.csv(glycolysis_stats, file.path(ct_dir, "glycolysis_genes_OA_vs_HC.csv"), row.names = FALSE)

  core_threshold_exp <- quantile(glycolysis_stats$mean_expr_oa, 0.75)
  core_genes <- glycolysis_stats[glycolysis_stats$mean_expr_oa >= core_threshold_exp &
                                   glycolysis_stats$logFC > 0, "gene"]

  message(ct, " - Core glycolysis genes (high expr in OA): ", length(core_genes))
  writeLines(as.character(core_genes), file.path(ct_dir, "core_glycolysis_genes.txt"))

  message("\n--- Virtual Knockout Analysis: ", ct, " ---")

  seu_ko <- seu_oa

  original_le_score <- seu_ko$LeadingEdgeLE_UCell
  original_le_ms <- seu_ko$LE_AddMS_1

  glycolysis_expr <- colMeans(expr_oa[intersect(ALL_GLYCOLYSIS_GENES, rownames(expr_oa)), ], na.rm = TRUE)
  glycolysis_score <- as.numeric(scale(glycolysis_expr))

  oa_scores <- data.frame(
    glycolysis = glycolysis_score,
    LE_UCell = original_le_score,
    LE_AddMS = original_le_ms
  )

  for (pathway in OA_PATHWAYS) {
    pathway_genes <- oa_pathway_genes[[pathway]]
    pathway_genes_present <- intersect(pathway_genes, rownames(expr_oa))
    if (length(pathway_genes_present) > 0) {
      pathway_score <- colMeans(expr_oa[pathway_genes_present, ], na.rm = TRUE)
      oa_scores[[pathway]] <- as.numeric(scale(pathway_score))
    }
  }

  write.csv(oa_scores, file.path(ct_dir, "original_scores.csv"), row.names = FALSE)

  ko_results <- list()
  for (gene in core_genes) {
    if (!(gene %in% rownames(expr_oa))) next

    expr_ko <- expr_oa
    expr_ko[gene, ] <- 0

    glycolysis_expr_ko <- colMeans(expr_ko[intersect(ALL_GLYCOLYSIS_GENES, rownames(expr_ko)), ], na.rm = TRUE)
    glycolysis_score_ko <- as.numeric(scale(glycolysis_expr_ko))

    ko_row <- data.frame(
      gene = gene,
      cell_type = ct,
      operation = "KO",
      glycolysis_score_orig = mean(glycolysis_score),
      glycolysis_score_ko = mean(glycolysis_score_ko),
      glycolysis_delta = mean(glycolysis_score_ko) - mean(glycolysis_score),
      LE_UCell_orig = mean(original_le_score),
      LE_UCell_ko = mean(original_le_score),
      LE_UCell_delta = 0,
      LE_AddMS_orig = mean(original_le_ms),
      LE_AddMS_ko = mean(original_le_ms),
      LE_AddMS_delta = 0
    )

    ko_results[[gene]] <- ko_row
  }

  ko_df <- dplyr::bind_rows(ko_results)
  write.csv(ko_df, file.path(ct_dir, "virtual_KO_results.csv"), row.names = FALSE)

  message("\n--- Virtual Overexpression Analysis: ", ct, " ---")

  oe_results <- list()
  for (gene in core_genes) {
    if (!(gene %in% rownames(expr_oa))) next

    expr_oe <- expr_oa
    expr_oe[gene, ] <- expr_oe[gene, ] * 2

    glycolysis_expr_oe <- colMeans(expr_oe[intersect(ALL_GLYCOLYSIS_GENES, rownames(expr_oe)), ], na.rm = TRUE)
    glycolysis_score_oe <- as.numeric(scale(glycolysis_expr_oe))

    oe_row <- data.frame(
      gene = gene,
      cell_type = ct,
      operation = "OE",
      glycolysis_score_orig = mean(glycolysis_score),
      glycolysis_score_oe = mean(glycolysis_score_oe),
      glycolysis_delta = mean(glycolysis_score_oe) - mean(glycolysis_score),
      LE_UCell_orig = mean(original_le_score),
      LE_UCell_oe = mean(original_le_score),
      LE_UCell_delta = 0,
      LE_AddMS_orig = mean(original_le_ms),
      LE_AddMS_oe = mean(original_le_ms),
      LE_AddMS_delta = 0
    )

    oe_results[[gene]] <- oe_row
  }

  oe_df <- dplyr::bind_rows(oe_results)
  write.csv(oe_df, file.path(ct_dir, "virtual_OE_results.csv"), row.names = FALSE)

  combined <- rbind(
    ko_df[, c("gene", "operation", "glycolysis_delta")],
    oe_df[, c("gene", "operation", "glycolysis_delta")]
  )
  combined$operation <- factor(combined$operation, levels = c("KO", "OE"))

  p_glycolysis <- ggplot2::ggplot(combined, ggplot2::aes(x = gene, y = glycolysis_delta, fill = operation)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge()) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste0("Virtual KO/OE Effects on Glycolysis - ", ct, " (OA)"),
      x = "Gene",
      y = "Change in Glycolysis Score",
      fill = "Operation"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))

  ggplot2::ggsave(file.path(ct_dir, "KO_OE_glycolysis_score_changes.png"), p_glycolysis, width = 10, height = 8)

  summary_stats <- data.frame(
    cell_type = ct,
    n_core_genes = length(core_genes),
    mean_glycolysis_score = mean(glycolysis_score),
    mean_LE_UCell = mean(original_le_score),
    mean_LE_AddMS = mean(original_le_ms),
    ko_mean_glycolysis_delta = mean(ko_df$glycolysis_delta),
    oe_mean_glycolysis_delta = mean(oe_df$glycolysis_delta),
    stringsAsFactors = FALSE
  )

  write.csv(summary_stats, file.path(ct_dir, "summary_statistics.csv"), row.names = FALSE)

  message("\n=== ", ct, " Analysis Complete ===")
  message("Core genes: ", length(core_genes))
  message("KO mean glycolysis delta: ", round(mean(ko_df$glycolysis_delta), 4))
  message("OE mean glycolysis delta: ", round(mean(oe_df$glycolysis_delta), 4))
}

message("\n========================================")
message("All analyses complete!")
message("Output directory: ", OUT_DIR)
message("========================================\n")

message("Summary of findings:")
for (ct in CELLTYPES) {
  ct_file <- file.path(OUT_DIR, ct, "summary_statistics.csv")
  if (file.exists(ct_file)) {
    summary <- read.csv(ct_file)
    message("\n", ct, ":")
    message("  - Core glycolysis genes: ", summary$n_core_genes)
    message("  - Mean glycolysis score: ", round(summary$mean_glycolysis_score, 3))
    message("  - KO effect on glycolysis: ", round(summary$ko_mean_glycolysis_delta, 4))
    message("  - OE effect on glycolysis: ", round(summary$oe_mean_glycolysis_delta, 4))
  }
}
