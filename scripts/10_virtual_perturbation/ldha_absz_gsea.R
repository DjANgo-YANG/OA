rm(list = ls())

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(fgsea)
  library(ggplot2)
})

ko_root <- Sys.getenv("OAGLY_KO_ROOT", unset = project_path("results_root", "virtual_ko"))
diff_file <- file.path(ko_root, "15_LDHA敲除", "diffRegulation_full.txt")
outdir <- Sys.getenv("OAGLY_LDHA_OUTDIR", unset = file.path(ko_root, "LDHA_HALLMARK_absZ"))
if (!file.exists(diff_file)) stop("LDHA KO input table was not found: ", diff_file)

diff <- fread(diff_file)
if (!all(c("gene", "Z") %in% names(diff))) stop("The LDHA KO table must contain gene and Z columns.")

# Retain the largest signed Z for duplicate symbols, as in the archived ranking.
rank_df <- diff %>%
  transmute(gene = toupper(trimws(gene)), rank_raw = as.numeric(Z)) %>%
  filter(!is.na(gene), gene != "", is.finite(rank_raw)) %>%
  group_by(gene) %>%
  summarise(rank_raw = max(rank_raw), .groups = "drop") %>%
  mutate(rank_abs = abs(rank_raw)) %>%
  arrange(desc(rank_abs), gene)
if (!nrow(rank_df) || !any(rank_df$rank_abs > 0)) stop("The LDHA KO table has no usable nonzero |Z| scores.")
stats <- setNames(rank_df$rank_abs, rank_df$gene)

gmt_file <- Sys.getenv("OAGLY_HALLMARK_GMT")
if (nzchar(gmt_file)) {
  pathways <- fgsea::gmtPathways(gmt_file)
} else {
  hallmark <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  pathways <- split(hallmark$gene_symbol, hallmark$gs_name)
}
pathways <- lapply(pathways, function(genes) unique(toupper(genes)))

set.seed(123)
enrichment <- fgsea::fgseaSimple(
  pathways = pathways,
  stats = stats,
  nperm = 10000,
  minSize = 10,
  maxSize = 500,
  scoreType = "pos",
  nproc = 1
) %>%
  as.data.frame() %>%
  arrange(pval, desc(NES), pathway) %>%
  mutate(leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";"))
if (!nrow(enrichment)) stop("No Hallmark pathways passed the size filters.")

top15 <- enrichment %>%
  slice_head(n = 15) %>%
  arrange(desc(NES))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
fwrite(rank_df, file.path(outdir, "GSEA_ranked_gene_by_absZ.csv"))
fwrite(enrichment, file.path(outdir, "HALLMARK_fgsea_absZ_results.csv"))
fwrite(filter(enrichment, padj < 0.05), file.path(outdir, "HALLMARK_fgsea_absZ_significant.csv"))
fwrite(top15, file.path(outdir, "LDHA_KO_HALLMARK_absZ_Top15.csv"))

plot_df <- top15 %>%
  mutate(label = gsub("_", " ", sub("^HALLMARK_", "", pathway)),
         label = factor(label, levels = rev(label)))
p <- ggplot(plot_df, aes(x = NES, y = label, fill = NES)) +
  geom_col(width = 0.72) +
  scale_fill_gradient(low = "#F2F8FD", high = "#245898", guide = "none") +
  labs(title = "LDHA KO | Hallmark enrichment by |Z|", x = "NES (|Z| enrichment)", y = NULL) +
  theme_classic()
ggsave(file.path(outdir, "LDHA_KO_HALLMARK_absZ_Top15.pdf"), p, width = 8, height = 6)

message("LDHA Hallmark enrichment by |Z| saved to: ", outdir)
