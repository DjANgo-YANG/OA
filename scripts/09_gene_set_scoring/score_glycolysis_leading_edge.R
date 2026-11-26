rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
options(stringsAsFactors = FALSE)

IN_RDS  <- project_path("results_root", "single_cell/annotation/seurat_object_celltype_0.6.rds")
OUT_DIR <- project_path("results_root", "glycolysis_scoring")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUT_DIR)

COL_DARK  <- "#08519c"
COL_LIGHT <- "#f7fbff"

pkgs <- c("Seurat", "dplyr", "ggplot2")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(Seurat)
library(dplyr)
library(ggplot2)

if (!requireNamespace("UCell", quietly = TRUE)) install.packages("UCell")
library(UCell)

obj <- readRDS(IN_RDS)
stopifnot("orig.ident" %in% colnames(obj@meta.data))
stopifnot("cell_type" %in% colnames(obj@meta.data))

obj$group <- ifelse(
  grepl("^OA", obj$orig.ident, ignore.case = TRUE), "OA",
  ifelse(grepl("^HC", obj$orig.ident, ignore.case = TRUE), "HC", NA)
)

cat("Group table:\n")
print(table(obj$group, useNA = "ifany"))
if (any(is.na(obj$group))) stop("group contains NA. Check orig.ident naming.")

leading_edge_str <- paste(
  "NT5E;MDH1;ANKZF1;CHST12;TALDO1;LDHA;CHST6;LCT;PKM;VEGFA;HSPA5;GNE;P4HA2;NDUFV3;PLOD2;CHPF2;SOD1;PFKFB1;AK3;PPP2CB;LHPP;ANGPTL4;B4GALT7;ENO2;FAM162A;ALDOA;SLC35A3;DPYSL4;AGL;CD44;AURKA;DCN;GOT1;IL13RA1;GAL3ST1;CHPF;NSDHL;QSOX1;HS2ST1;ECD;TPST1;SPAG4;ERO1A;IGFBP3;PGM2;ELF3;B4GALT1;PPFIA4;GPC1;NANP;GPC3;FBP2;NASP;SDHC;KDELR3;TPI1;STC1;XYLT2;B4GALT2;SRD5A3;ENO1;B4GALT4;EXT1;MED24;GMPPB;PPIA;AK4;GALE;PC;CASP6;GYS1;CHST4;GFUS;PYGL;PAM;SLC37A4;PYGB;MPI;GLCE;MERTK",
  sep = ""
)
gene_set_raw <- unique(trimws(unlist(strsplit(leading_edge_str, ";"))))
gene_set_raw <- gene_set_raw[gene_set_raw != ""]

feat <- rownames(obj)

to_title <- function(x) paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
gene_set_try <- list(
  raw   = gene_set_raw,
  upper = toupper(gene_set_raw),
  title = to_title(tolower(gene_set_raw))
)

match_counts <- sapply(gene_set_try, function(gs) sum(gs %in% feat))
best_key <- names(which.max(match_counts))
gene_set <- gene_set_try[[best_key]]

genes_in  <- gene_set[gene_set %in% feat]
genes_out <- setdiff(gene_set_raw, genes_in)

cat("Gene matching mode:", best_key, "\n")
cat("Provided:", length(gene_set_raw), " Found:", length(genes_in), " Missing:", length(genes_out), "\n")

write.table(
  data.frame(gene = gene_set_raw, in_object = gene_set_raw %in% genes_in),
  file = "LeadingEdge_gene_presence.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

if (length(genes_in) < 10) stop("Too few genes matched (<10). Check gene IDs (Symbol vs Ensembl).")

obj <- AddModuleScore(obj, features = list(LeadingEdge = genes_in), name = "LE_AddMS_")
addms_col <- grep("^LE_AddMS_", colnames(obj@meta.data), value = TRUE)[1]
cat("Using AddModuleScore column:", addms_col, "\n")

obj <- UCell::AddModuleScore_UCell(obj, features = list(LeadingEdge = genes_in), name = "LE_UCell")

ucell_col <- grep("UCell$", colnames(obj@meta.data), value = TRUE)
if (length(ucell_col) == 0) ucell_col <- grep("UCell", colnames(obj@meta.data), value = TRUE)
if (length(ucell_col) == 0) stop("No UCell score column found. UCell scoring may have failed.")

if (length(ucell_col) > 1) {
  cat("Multiple UCell columns found, using the first one:\n")
  print(ucell_col)
}
ucell_col <- ucell_col[1]
cat("Using UCell score column:", ucell_col, "\n")

plot_global <- function(score_col, label) {
  stopifnot(score_col %in% colnames(obj@meta.data))

  df <- obj@meta.data %>%
    dplyr::select(group, cell_type, all_of(score_col)) %>%
    dplyr::rename(score = all_of(score_col)) %>%
    dplyr::filter(group %in% c("HC", "OA"))

  pval <- wilcox.test(score ~ group, data = df)$p.value

  g <- ggplot(df, aes(x = group, y = score)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.2, outlier.size = 0.2) +
    labs(
      title = paste0(label, " (global)"),
      subtitle = paste0("Wilcoxon p = ", signif(pval, 3)),
      x = NULL,
      y = "Score"
    ) +
    theme_classic(base_size = 12)

  ggsave(paste0(label, "_Global_violin.png"), g, width = 5.5, height = 4.2, dpi = 300)
  data.frame(method = label, level = "global", p_value = pval)
}

plot_by_celltype <- function(score_col, label) {
  stopifnot(score_col %in% colnames(obj@meta.data))

  df <- obj@meta.data %>%
    dplyr::select(group, cell_type, all_of(score_col)) %>%
    dplyr::rename(score = all_of(score_col)) %>%
    dplyr::filter(group %in% c("HC", "OA"))

  celltypes <- sort(unique(df$cell_type))

  res <- lapply(celltypes, function(ct) {
    d <- df %>% dplyr::filter(cell_type == ct)
    nHC <- sum(d$group == "HC")
    nOA <- sum(d$group == "OA")

    if (nHC < 20 || nOA < 20) {
      return(data.frame(
        cell_type = ct,
        n_HC = nHC,
        n_OA = nOA,
        mean_HC = mean(d$score[d$group == "HC"], na.rm = TRUE),
        mean_OA = mean(d$score[d$group == "OA"], na.rm = TRUE),
        delta = NA_real_,
        p_value = NA_real_
      ))
    }

    p <- wilcox.test(score ~ group, data = d)$p.value
    mHC <- mean(d$score[d$group == "HC"], na.rm = TRUE)
    mOA <- mean(d$score[d$group == "OA"], na.rm = TRUE)

    g <- ggplot(d, aes(x = group, y = score)) +
      geom_violin(trim = TRUE) +
      geom_boxplot(width = 0.2, outlier.size = 0.2) +
      labs(
        title = paste0(label, " | ", ct),
        subtitle = paste0("n(HC)=", nHC, "  n(OA)=", nOA, "  p=", signif(p, 3)),
        x = NULL,
        y = "Score"
      ) +
      theme_classic(base_size = 11)

    safe_ct <- gsub("[^A-Za-z0-9_\\-]", "_", ct)
    ggsave(paste0(label, "_", safe_ct, ".png"), g, width = 5.5, height = 4.2, dpi = 300)

    data.frame(
      cell_type = ct,
      n_HC = nHC,
      n_OA = nOA,
      mean_HC = mHC,
      mean_OA = mOA,
      delta = mOA - mHC,
      p_value = p
    )
  }) %>% bind_rows()

  res$method <- label
  res$p_adj <- p.adjust(res$p_value, method = "BH")

  write.table(
    res,
    paste0(label, "_ByCellType_stats.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  gsum_df <- res %>%
    dplyr::filter(!is.na(delta)) %>%
    dplyr::mutate(
      cell_type = reorder(cell_type, delta),
      fill_value = rank(abs(delta), ties.method = "first")
    )

  gsum <- ggplot(gsum_df, aes(x = cell_type, y = delta, fill = fill_value)) +
    geom_col(width = 0.8) +
    coord_flip() +
    scale_fill_gradient(low = COL_LIGHT, high = COL_DARK) +
    labs(
      title = paste0(label, " Δ(OA-HC) by cell type"),
      x = NULL,
      y = "Mean score difference (OA - HC)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )

  ggsave(
    paste0(label, "_DeltaBar_byCellType.png"),
    gsum,
    width = 7.5,
    height = 5.5,
    dpi = 300
  )

  res
}

stats_global <- bind_rows(
  plot_global(addms_col, "LeadingEdge_AddModuleScore"),
  plot_global(ucell_col, "LeadingEdge_UCell")
)

res_add <- plot_by_celltype(addms_col, "LeadingEdge_AddModuleScore")
res_uc  <- plot_by_celltype(ucell_col, "LeadingEdge_UCell")

write.table(
  stats_global,
  "LeadingEdge_global_stats.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

saveRDS(obj, file = "seurat_with_LeadingEdge_scores.rds")

cat("\nDONE. Output dir:", OUT_DIR, "\n")
