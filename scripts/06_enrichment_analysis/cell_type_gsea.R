rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
options(stringsAsFactors = FALSE)

INPUT_RDS <- project_path("results_root", "single_cell/annotation/seurat_object_celltype_0.6.rds")
OUT_DIR   <- project_path("results_root", "cell_type_enrichment")

SAMPLE_COL   <- "orig.ident"
CELLTYPE_COL <- "cell_type"
NEW_GROUP    <- "condition"

ASSAY_USE <- "RNA"
SLOT_USE  <- "counts"

MIN_CELLS_PER_SAMPLE  <- 30
MIN_SAMPLES_PER_GROUP <- 3

SPECIES <- "Homo sapiens"
P_ADJ_CUTOFF <- 0.05
SHOW_TOP_N <- 20
KEYWORDS <- c("EPITHELIAL_MESENCHYMAL","TNFA","HYPOXIA","MATRIX","COLLAGEN","ECM",
              "INFLAM","IL","WNT","TGF","ANGIO","APOPT","GLYCOL")

set.seed(1234)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

pkgs_cran <- c("ggplot2", "dplyr", "tibble")
pkgs_bioc <- c("Seurat", "msigdbr", "clusterProfiler", "enrichplot", "limma", "edgeR")

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

if (!SAMPLE_COL %in% colnames(seu@meta.data)) stop(paste0("Cannot find SAMPLE_COL: ", SAMPLE_COL))
if (!CELLTYPE_COL %in% colnames(seu@meta.data)) stop(paste0("Cannot find CELLTYPE_COL: ", CELLTYPE_COL))
if (!ASSAY_USE %in% names(seu@assays)) stop(paste0("Cannot find assay: ", ASSAY_USE))

DefaultAssay(seu) <- ASSAY_USE

hc_samples <- paste0("HC", 1:6)
oa_samples <- paste0("OA", 1:6)

all_samples <- unique(as.character(seu@meta.data[[SAMPLE_COL]]))
missing_hc <- setdiff(hc_samples, all_samples)
missing_oa <- setdiff(oa_samples, all_samples)
if (length(missing_hc) > 0 || length(missing_oa) > 0) {
  stop(paste0(
    "orig.ident missing expected sample names.\n",
    "missing HC: ", paste(missing_hc, collapse = ", "), "\n",
    "missing OA: ", paste(missing_oa, collapse = ", "), "\n",
    "Actual orig.ident: ", paste(sort(all_samples), collapse = ", ")
  ))
}

seu@meta.data[[NEW_GROUP]] <- NA_character_
seu@meta.data[[NEW_GROUP]][seu@meta.data[[SAMPLE_COL]] %in% hc_samples] <- "HC"
seu@meta.data[[NEW_GROUP]][seu@meta.data[[SAMPLE_COL]] %in% oa_samples] <- "OA"

if (any(is.na(seu@meta.data[[NEW_GROUP]]))) {
  bad <- unique(as.character(seu@meta.data[[SAMPLE_COL]][is.na(seu@meta.data[[NEW_GROUP]])]))
  stop(paste0("Unmapped orig.ident values: ", paste(bad, collapse = ", "),
              ". Fix mapping or filter them out."))
}

message("Cells by condition:")
print(table(seu@meta.data[[NEW_GROUP]]))
message("Cells by sample:")
print(table(seu@meta.data[[SAMPLE_COL]]))

message("Preparing MSigDB gene sets (Hallmark + GO:BP) ...")
m_H  <- msigdbr::msigdbr(species = SPECIES, category = "H") %>% dplyr::select(gs_name, gene_symbol)
m_C5 <- msigdbr::msigdbr(species = SPECIES, category = "C5", subcategory = "GO:BP") %>% dplyr::select(gs_name, gene_symbol)
m_t2g <- rbind(m_H, m_C5)

celltypes <- sort(unique(as.character(seu@meta.data[[CELLTYPE_COL]])))
message("Found cell types: ", length(celltypes))

run_log <- list()

for (ct in celltypes) {
  message("\n==============================")
  message("Cell type: ", ct)
  message("==============================")

  ct_dir <- file.path(OUT_DIR, gsub("[^A-Za-z0-9_\\-]+", "_", ct))
  if (!dir.exists(ct_dir)) dir.create(ct_dir, recursive = TRUE)

  cells_ct <- rownames(seu@meta.data)[seu@meta.data[[CELLTYPE_COL]] == ct]
  if (length(cells_ct) < 50) {
    message("Skip: too few cells overall (", length(cells_ct), ")")
    run_log[[ct]] <- data.frame(cell_type=ct, status="skip_few_cells_total", stringsAsFactors = FALSE)
    next
  }
  seu_ct <- subset(seu, cells = cells_ct)

  meta_ct <- seu_ct@meta.data
  tab_ct <- sort(table(meta_ct[[SAMPLE_COL]]), decreasing = TRUE)
  write.csv(as.data.frame(tab_ct), file.path(ct_dir, "cells_per_sample_in_celltype.csv"))

  keep_samples <- names(tab_ct)[tab_ct >= MIN_CELLS_PER_SAMPLE]
  keep_hc <- intersect(keep_samples, hc_samples)
  keep_oa <- intersect(keep_samples, oa_samples)

  if (length(keep_hc) < MIN_SAMPLES_PER_GROUP || length(keep_oa) < MIN_SAMPLES_PER_GROUP) {
    message("Skip: group samples insufficient after filtering. HC=", length(keep_hc), ", OA=", length(keep_oa))
    run_log[[ct]] <- data.frame(cell_type=ct, status="skip_group_unbalanced_after_filter", stringsAsFactors = FALSE)
    next
  }

  cells_keep <- rownames(seu_ct@meta.data)[seu_ct@meta.data[[SAMPLE_COL]] %in% c(keep_hc, keep_oa)]
  seu_ct <- subset(seu_ct, cells = cells_keep)

  message("Pseudobulk aggregation ...")
  pb_list <- AggregateExpression(
    seu_ct,
    group.by = SAMPLE_COL,
    assays = ASSAY_USE,
    slot = SLOT_USE,
    verbose = FALSE
  )
  pb_counts <- round(pb_list[[ASSAY_USE]])
  pb_counts <- pb_counts[, colnames(pb_counts) %in% c(keep_hc, keep_oa), drop = FALSE]

  pb_counts <- pb_counts[, c(sort(intersect(colnames(pb_counts), keep_hc)),
                             sort(intersect(colnames(pb_counts), keep_oa))), drop = FALSE]

  write.csv(as.data.frame(pb_counts), file.path(ct_dir, "pseudobulk_counts_genes_by_sample.csv"))

  sample_info <- data.frame(
    sample = colnames(pb_counts),
    condition = ifelse(colnames(pb_counts) %in% keep_oa, "OA", "HC"),
    stringsAsFactors = FALSE
  )
  sample_info$condition <- factor(sample_info$condition, levels = c("HC", "OA"))
  write.csv(sample_info, file.path(ct_dir, "sample_info.csv"), row.names = FALSE)

  message("DE (limma-voom) ...")
  dge <- edgeR::DGEList(counts = pb_counts, samples = sample_info)
  keep_gene <- edgeR::filterByExpr(dge, group = sample_info$condition)
  dge <- dge[keep_gene, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)

  design <- model.matrix(~ 0 + condition, data = dge$samples)
  colnames(design) <- gsub("^condition", "", colnames(design))
  contrast <- limma::makeContrasts(OA_vs_HC = OA - HC, levels = design)

  v <- limma::voom(dge, design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  fit2 <- limma::contrasts.fit(fit, contrast)
  fit2 <- limma::eBayes(fit2)

  de_tab <- limma::topTable(fit2, coef = "OA_vs_HC", number = Inf, sort.by = "none")
  de_tab$gene <- rownames(de_tab)
  de_tab <- de_tab %>% dplyr::select(gene, logFC, AveExpr, t, P.Value, adj.P.Val, B)

  write.csv(de_tab, file.path(ct_dir, "DE_limma_voom_OA_vs_HC_full.csv"), row.names = FALSE)
  de_sig <- de_tab %>% dplyr::filter(adj.P.Val < P_ADJ_CUTOFF)
  write.csv(de_sig, file.path(ct_dir, "DE_limma_voom_OA_vs_HC_sig.csv"), row.names = FALSE)

  gene_list <- de_tab$t
  names(gene_list) <- de_tab$gene
  gene_list <- sort(gene_list, decreasing = TRUE)
  write.csv(
    data.frame(Gene = names(gene_list), t_stat = as.numeric(gene_list)),
    file.path(ct_dir, "Ranked_Gene_List_tstat.csv"),
    row.names = FALSE
  )

  message("GSEA ...")
  gsea_res <- clusterProfiler::GSEA(
    geneList     = gene_list,
    TERM2GENE    = m_t2g,
    pvalueCutoff = 1,
    minGSSize    = 10,
    maxGSSize    = 500,
    eps          = 0,
    seed         = TRUE
  )

  gsea_df <- gsea_res@result
  write.csv(gsea_df, file.path(ct_dir, "GSEA_full_results.csv"), row.names = FALSE)

  gsea_sig <- gsea_df %>% dplyr::filter(p.adjust < P_ADJ_CUTOFF)
  write.csv(gsea_sig, file.path(ct_dir, "GSEA_sig_FDRlt0.05.csv"), row.names = FALSE)

  if (nrow(gsea_sig) > 0) {
    message("Plots ...")
    p_ridge <- enrichplot::ridgeplot(gsea_res, showCategory = min(SHOW_TOP_N, nrow(gsea_sig))) +
      ggplot2::labs(title = paste0("GSEA Ridgeplot (", ct, ") OA vs HC, pseudobulk rank=t"))
    ggplot2::ggsave(file.path(ct_dir, "GSEA_ridgeplot.png"), p_ridge, width = 10, height = 12)

    top_up   <- gsea_sig %>% dplyr::arrange(dplyr::desc(NES)) %>% head(10) %>% dplyr::pull(ID)
    top_down <- gsea_sig %>% dplyr::arrange(NES) %>% head(10) %>% dplyr::pull(ID)
    top_paths <- unique(c(top_up, top_down))

    p_dot <- enrichplot::dotplot(gsea_res, showCategory = top_paths, split = ".sign") +
      ggplot2::facet_grid(. ~ .sign) +
      ggplot2::ggtitle(paste0("Top Pathways (", ct, ") OA vs HC - pseudobulk")) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
    ggplot2::ggsave(file.path(ct_dir, "GSEA_top_dotplot.png"), p_dot, width = 12, height = 10)

    hit_ids <- grep(paste(KEYWORDS, collapse="|"), gsea_sig$ID, value=TRUE, ignore.case=TRUE)
    if (length(hit_ids) > 0) {
      for (pid in head(hit_ids, 10)) {
        p <- enrichplot::gseaplot2(gsea_res, geneSetID = pid, title = pid, pvalue_table = TRUE)
        fname <- gsub("[:/\\\\ ]", "_", pid)
        ggplot2::ggsave(file.path(ct_dir, paste0("GSEA_", fname, ".png")), p, width = 8, height = 6)
      }
    }
  }

  run_log[[ct]] <- data.frame(
    cell_type = ct,
    status = "ok",
    n_cells = ncol(seu_ct),
    n_samples_used = ncol(pb_counts),
    n_HC = sum(sample_info$condition == "HC"),
    n_OA = sum(sample_info$condition == "OA"),
    n_sig_genes = nrow(de_sig),
    n_sig_paths = nrow(gsea_sig),
    stringsAsFactors = FALSE
  )

  message("Done: ", ct, " | sig genes=", nrow(de_sig), " | sig pathways=", nrow(gsea_sig))
}

log_df <- dplyr::bind_rows(run_log)
write.csv(log_df, file.path(OUT_DIR, "RUN_SUMMARY_by_celltype.csv"), row.names = FALSE)

message("\nAll done. Output: ", OUT_DIR)
message("Key summary file: RUN_SUMMARY_by_celltype.csv")
