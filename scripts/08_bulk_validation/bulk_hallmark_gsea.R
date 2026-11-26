source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))

options(stringsAsFactors = FALSE)

WORKDIR    <- project_path("data_root", "bulk")
COUNT_FILE <- file.path(WORKDIR, "GSE254844_processed_counts.txt")
OUTDIR     <- file.path(WORKDIR, "Bulk_OA_vs_HC_UnbiasedGSEA")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
setwd(OUTDIR)

EXP_FRAC <- 0.20

MIN_SIZE <- 10
MAX_SIZE <- 500
SEED     <- 1

SPECIES  <- "Homo sapiens"

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(dplyr)
  library(tibble)
  library(msigdbr)
  library(fgsea)
  library(ggplot2)
})

cts <- read.table(
  COUNT_FILE,
  header       = TRUE,
  sep          = "",
  quote        = "",
  comment.char = "",
  check.names  = FALSE,
  row.names    = 1
)
write.csv(cts, "QC_counts_raw.csv")

samples <- colnames(cts)

condition <- ifelse(grepl("^Control\\.", samples), "HC",
                    ifelse(grepl("^(WB|NWB)\\.", samples), "OA", NA))
if (any(is.na(condition))) {
  bad <- samples[is.na(condition)]
  stop("Condition parsing failed. Example bad samples: ", paste(head(bad, 10), collapse = ", "))
}

meta <- data.frame(sample_id = samples, condition = factor(condition, levels = c("HC","OA")))
write.csv(meta, "sample_info_OA_vs_HC.csv", row.names = FALSE)
cat("Sample counts:\n")
print(table(meta$condition))

keep <- rowMeans(cts > 0) >= EXP_FRAC
cts_f <- cts[keep, , drop = FALSE]
cat("Genes kept after expression filter:", nrow(cts_f), "\n")

dge <- DGEList(counts = cts_f)
dge <- calcNormFactors(dge)

design <- model.matrix(~ 0 + condition, data = meta)
colnames(design) <- c("HC","OA")

v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- contrasts.fit(fit, makeContrasts(OA_vs_HC = OA - HC, levels = design))
fit <- eBayes(fit, trend = TRUE)

tt <- topTable(fit, coef = "OA_vs_HC", number = Inf, sort.by = "none") %>%
  rownames_to_column("gene")
write.csv(tt, "DE_OA_vs_HC.csv", row.names = FALSE)

ranks <- tt$t
names(ranks) <- tt$gene
ranks <- ranks[!is.na(ranks)]
ranks <- sort(ranks, decreasing = TRUE)

genes_universe <- names(ranks)

cols <- msigdbr_collections()
write.csv(cols, "msigdbr_collections_available.csv", row.names = FALSE)

get_pathways_safe <- function(species, category, subcategory = NULL,
                              genes_universe, min_size = 10) {
  cols0 <- msigdbr_collections()

  if (!is.null(subcategory)) {
    ok <- any(cols0$category == category & cols0$subcategory == subcategory)
    if (!ok) {
      message("Skip: ", category, " / ", subcategory, " (NOT available in your msigdbr).")
      return(NULL)
    }
    ms <- msigdbr(species = species, category = category, subcategory = subcategory)
  } else {

    ok <- any(cols0$category == category) || any(cols0$collection == category)

    ms <- tryCatch(
      msigdbr(species = species, collection = category),
      error = function(e) NULL
    )
    if (is.null(ms)) {
      message("Skip: collection=", category, " (NOT available in your msigdbr).")
      return(NULL)
    }
  }

  ms <- ms %>% dplyr::select(gs_name, gene_symbol)
  pw <- split(ms$gene_symbol, ms$gs_name)
  pw <- lapply(pw, function(gs) intersect(gs, genes_universe))
  pw <- pw[lengths(pw) >= min_size]
  pw
}

pathways_H    <- get_pathways_safe(SPECIES, category = "H", genes_universe = genes_universe, min_size = MIN_SIZE)
pathways_RE   <- get_pathways_safe(SPECIES, category = "C2", subcategory = "CP:REACTOME", genes_universe = genes_universe, min_size = MIN_SIZE)
pathways_KEGG <- get_pathways_safe(SPECIES, category = "C2", subcategory = "CP:KEGG",     genes_universe = genes_universe, min_size = MIN_SIZE)
pathways_GOBP <- get_pathways_safe(SPECIES, category = "C5", subcategory = "GO:BP",       genes_universe = genes_universe, min_size = MIN_SIZE)

path_list <- list(
  Hallmark = pathways_H,
  Reactome = pathways_RE,
  KEGG     = pathways_KEGG,
  GOBP     = pathways_GOBP
)
path_list <- path_list[!vapply(path_list, is.null, logical(1))]

cat("\nPathway collections to run:\n")
print(sapply(path_list, length))

if (length(path_list) == 0) {
  stop("No pathway collections available. Check msigdbr_collections_available.csv and your SPECIES setting.")
}

set.seed(SEED)

run_fgsea <- function(pathways, ranks, label, out_csv) {
  fg <- fgseaMultilevel(pathways = pathways, stats = ranks,
                        minSize = MIN_SIZE, maxSize = MAX_SIZE)
  fg <- fg[order(fg$padj, -abs(fg$NES)), ]
  fg_df <- as.data.frame(fg)
  fg_df$collection <- label
  if ("leadingEdge" %in% colnames(fg_df)) {
    fg_df$leadingEdge <- vapply(fg_df$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  }
  write.csv(fg_df, out_csv, row.names = FALSE)
  fg_df
}

fg_all_list <- list()
for (nm in names(path_list)) {
  out_csv <- paste0("GSEA_", nm, "_all.csv")
  fg_all_list[[nm]] <- run_fgsea(path_list[[nm]], ranks, nm, out_csv)
}

fg_all <- bind_rows(fg_all_list)
write.csv(fg_all, "GSEA_ALL_collections.csv", row.names = FALSE)

met_kw <- c(
  "METABOL", "GLYCOL", "GLUC", "PYRUV", "LACT", "TCA", "CITRIC",
  "OXIDAT", "RESPIR", "ELECTRON", "MITOCHON", "ATP", "NADH",
  "FATTY", "LIPID", "BETA_OXID", "CHOLESTER", "PENTOSE",
  "AMINO", "NUCLEOT", "PURINE", "PYRIMID", "CARBOHYDR",
  "HEXOS", "FRUCT", "GALACT", "PHOSPHATE", "REDOX"
)

met_all <- fg_all %>%
  dplyr::filter(grepl(paste(met_kw, collapse = "|"), pathway, ignore.case = TRUE)) %>%
  dplyr::arrange(padj, desc(abs(NES)))

write.csv(met_all, "GSEA_metabolism_related_posthoc.csv", row.names = FALSE)
write.csv(met_all %>% dplyr::slice_head(n = 50),
          "GSEA_metabolism_related_TOP50.csv", row.names = FALSE)

plot_top <- function(df, title, out_png, n = 25) {
  df2 <- df %>%
    dplyr::filter(!is.na(padj), !is.na(NES)) %>%
    dplyr::slice_head(n = n) %>%
    dplyr::mutate(pathway = factor(pathway, levels = rev(pathway)))

  p <- ggplot(df2, aes(x = pathway, y = NES, fill = NES)) +
    geom_col(width = 0.8) +
    coord_flip() +
    scale_fill_gradient(
      low = "#f7fbff",
      high = "#08519c"
    ) +
    theme_bw() +
    labs(title = title, x = NULL, y = "NES (OA vs HC)") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )

  ggsave(out_png, p, width = 10, height = 7)
}

plot_top(fg_all %>% dplyr::arrange(padj, desc(abs(NES))),
         "Unbiased GSEA (ALL collections): top NES/padj",
         "GSEA_ALL_top25.png", n = 25)

plot_top(met_all,
         "Unbiased GSEA (post hoc metabolism-related): top NES/padj",
         "Metabolism_posthoc_top25.png", n = 25)

cat("\nTop 10 overall pathways (by padj, |NES|):\n")
print(fg_all %>% dplyr::select(collection, pathway, NES, padj, size) %>% dplyr::slice_head(n = 10))

cat("\nTop 10 metabolism-related pathways (post hoc):\n")
print(met_all %>% dplyr::select(collection, pathway, NES, padj, size) %>% dplyr::slice_head(n = 10))

cat("\nDONE.\nOutputs in: ", OUTDIR, "\n",
    "Key files:\n",
    " - msigdbr_collections_available.csv\n",
    " - DE_OA_vs_HC.csv\n",
    " - GSEA_ALL_collections.csv\n",
    " - GSEA_metabolism_related_posthoc.csv (+TOP50)\n",
    " - GSEA_ALL_top25.png\n",
    " - Metabolism_posthoc_top25.png\n", sep = "")
