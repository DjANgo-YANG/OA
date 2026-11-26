rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
options(stringsAsFactors = FALSE)

cran_pkgs <- c("data.table", "dplyr", "stringr", "tidyr", "pheatmap")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(pheatmap)

base_dir <- project_path("results_root", "virtual_oe")
gsea_dir <- file.path(base_dir, "HALLMARK_GSEA_EnrichmentMap")
out_dir  <- file.path(base_dir, "HALLMARK_GSEA_heatmap")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

all_gsea <- fread(file.path(gsea_dir, "All_OE_HALLMARK_GSEA_all.csv"))

heat_input <- all_gsea %>%
  transmute(
    folder = paste0(sprintf("%02d", match(OE_gene, unique(OE_gene))), "_", OE_gene, "_", OE_gene, "过表达"),
    OE     = OE_gene,
    gs_name = Description,
    NES    = NES,
    padj   = FDR
  )

fwrite(heat_input, file.path(out_dir, "HALLMARK_OE_heatmap_input.csv"))

nes_mat <- heat_input %>%
  select(OE, gs_name, NES) %>%
  pivot_wider(names_from = OE, values_from = NES, values_fn = mean, values_fill = 0) %>%
  as.data.frame()

padj_mat <- heat_input %>%
  select(OE, gs_name, padj) %>%
  pivot_wider(names_from = OE, values_from = padj, values_fn = min, values_fill = 1) %>%
  as.data.frame()

rownames(nes_mat)  <- nes_mat$gs_name
rownames(padj_mat) <- padj_mat$gs_name
nes_mat$gs_name  <- NULL
padj_mat$gs_name <- NULL

clean_names <- rownames(nes_mat) %>%
  str_replace_all("^HALLMARK_", "") %>%
  str_replace_all("_", " ")
rownames(nes_mat)  <- clean_names
rownames(padj_mat) <- clean_names

oe_genes <- all_gsea %>% pull(OE_gene) %>% unique()
oe_genes <- oe_genes[order(match(oe_genes, unique(all_gsea$OE_gene)))]
nes_mat  <- nes_mat[, oe_genes, drop = FALSE]
padj_mat <- padj_mat[, oe_genes, drop = FALSE]

sig_mat <- padj_mat
sig_mat[padj_mat < 0.05]  <- "*"
sig_mat[padj_mat >= 0.05] <- ""

palette_len <- 100
my_palette <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(palette_len)

pdf(file.path(out_dir, "HALLMARK_23OE_heatmap_NES.pdf"), width = 14, height = 10)
pheatmap(
  nes_mat,
  color = my_palette,
  breaks = seq(-max(abs(nes_mat), na.rm = TRUE), max(abs(nes_mat), na.rm = TRUE), length.out = palette_len + 1),
  display_numbers = sig_mat,
  fontsize_number = 8,
  number_color = "black",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 9,
  fontsize_col = 9,
  angle_col = 45,
  main = "OE Hallmark GSEA heatmap (23 OE genes, NES)",
  legend_labels = "NES",
  border_color = "grey90"
)
dev.off()

png(file.path(out_dir, "HALLMARK_23OE_heatmap_NES.png"), width = 14, height = 10,
    units = "in", res = 300)
pheatmap(
  nes_mat,
  color = my_palette,
  breaks = seq(-max(abs(nes_mat), na.rm = TRUE), max(abs(nes_mat), na.rm = TRUE), length.out = palette_len + 1),
  display_numbers = sig_mat,
  fontsize_number = 8,
  number_color = "black",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 9,
  fontsize_col = 9,
  angle_col = 45,
  main = "OE Hallmark GSEA heatmap (23 OE genes, NES)",
  legend_labels = "NES",
  border_color = "grey90"
)
dev.off()

message("完成。输出目录: ", out_dir)
