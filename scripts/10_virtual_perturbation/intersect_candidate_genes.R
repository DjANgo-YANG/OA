rm(list = ls())
gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
options(stringsAsFactors = FALSE)

cran_pkgs <- c("data.table", "dplyr", "stringr", "tibble")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c("Seurat", "edgeR", "limma")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
}

library(data.table)
library(dplyr)
library(stringr)
library(tibble)
library(Seurat)
library(edgeR)
library(limma)
if (!requireNamespace("ggalluvial", quietly = TRUE)) install.packages("ggalluvial")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")

library(ggalluvial)
library(ggplot2)
library(RColorBrewer)
if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
library(VennDiagram)
library(grid)

ko_dirs <- c(
  project_path("data_root", "figure_inputs/virtual_ko/20_PPIA敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/21_AK3敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/22_TGFBI敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/23_COG2敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/01_AGL敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/02_NASP敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/03_LHX9敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/04_COPB2敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/05_PAM敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/06_VCAN敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/07_ELF3敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/08_SDC2敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/09_ME2敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/10_PGK1敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/11_HDLBP敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/12_MDH1敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/13_DCN敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/14_PKM敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/15_LDHA敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/16_SLC35A3敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/17_PLOD2敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/18_CD44敲除"),
  project_path("data_root", "figure_inputs/virtual_ko/19_ENO1敲除")
)

seu_path  <- project_path("results_root", "single_cell/annotation/seurat_object_celltype_0.6.rds")
bulk_path <- project_path("data_root", "bulk/OA_COUNTS_raw.csv")

out_dir <- project_path("data_root", "figure_inputs/virtual_ko/KO_intersection_with_sc_and_bulk")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sc_logfc_cutoff <- 0.25
sc_padj_cutoff  <- 0.05
sc_min_pct      <- 0.10

bulk_logfc_cutoff <- 0.5
bulk_padj_cutoff  <- 0.05

extract_ko_gene <- function(x) {
  x <- basename(x)
  x <- str_remove(x, "^\\d+_")
  x <- str_remove(x, "敲除$")
  x
}

ko_genes <- unique(extract_ko_gene(ko_dirs))
ko_df <- data.frame(KO_gene = ko_genes, stringsAsFactors = FALSE)

write.csv(ko_df, file.path(out_dir, "KO_genes_extracted.csv"), row.names = FALSE)

cat("提取到 KO genes:\n")
print(ko_genes)

seu <- readRDS(seu_path)

table(seu@meta.data$orig.ident)

seu@meta.data$group <- ifelse(
  grepl("^HC", seu@meta.data$orig.ident),
  "HC",
  ifelse(grepl("^OA", seu@meta.data$orig.ident), "OA", NA)
)

seu@meta.data$group <- factor(seu@meta.data$group, levels = c("HC", "OA"))

table(seu@meta.data$group, useNA = "ifany")
table(seu@meta.data$orig.ident, seu@meta.data$group, useNA = "ifany")
if (!"cell_type" %in% colnames(seu@meta.data)) {
  stop("Seurat对象meta.data中未找到 cell_type 列")
}
if (!"group" %in% colnames(seu@meta.data)) {
  stop("Seurat对象meta.data中未找到 group 列")
}

cat("\n单细胞 group 分布：\n")
print(table(seu@meta.data$group, useNA = "ifany"))

cat("\n单细胞 cell_type 分布：\n")
print(table(seu@meta.data$cell_type, useNA = "ifany"))

Idents(seu) <- seu@meta.data$group

run_sc_deg_one_celltype <- function(seu_obj, ct_name,
                                    group_col = "group",
                                    logfc_cutoff = 0.25,
                                    padj_cutoff = 0.05,
                                    min_pct = 0.1) {
  sub <- subset(seu_obj, subset = cell_type == ct_name)

  if (ncol(sub) == 0) {
    warning("cell_type ", ct_name, " 无细胞，跳过")
    return(NULL)
  }

  grp_tab <- table(sub@meta.data[[group_col]])
  cat("\n", ct_name, " group分布:\n", sep = "")
  print(grp_tab)

  if (!all(c("HC", "OA") %in% names(grp_tab))) {
    warning(ct_name, " 中缺少 HC 或 OA，跳过")
    return(NULL)
  }

  Idents(sub) <- sub@meta.data[[group_col]]

  deg <- FindMarkers(
    object = sub,
    ident.1 = "OA",
    ident.2 = "HC",
    logfc.threshold = 0,
    min.pct = min_pct,
    test.use = "wilcox"
  )

  if (nrow(deg) == 0) return(NULL)

  deg <- deg %>%
    rownames_to_column("gene")

  logfc_col <- if ("avg_log2FC" %in% colnames(deg)) "avg_log2FC" else
    if ("avg_logFC" %in% colnames(deg)) "avg_logFC" else
      stop("FindMarkers结果中未找到 avg_log2FC / avg_logFC 列")

  deg <- deg %>%
    mutate(
      cell_type = ct_name,
      logFC_use = .data[[logfc_col]]
    ) %>%
    filter(!is.na(p_val_adj), p_val_adj < padj_cutoff, abs(logFC_use) >= logfc_cutoff)

  deg
}

deg_preHTC <- run_sc_deg_one_celltype(
  seu_obj = seu,
  ct_name = "preHTC",
  logfc_cutoff = sc_logfc_cutoff,
  padj_cutoff = sc_padj_cutoff,
  min_pct = sc_min_pct
)

deg_HTC <- run_sc_deg_one_celltype(
  seu_obj = seu,
  ct_name = "HTC",
  logfc_cutoff = sc_logfc_cutoff,
  padj_cutoff = sc_padj_cutoff,
  min_pct = sc_min_pct
)

deg_sc_all <- bind_rows(deg_preHTC, deg_HTC)

if (nrow(deg_sc_all) == 0) {
  warning("单细胞 preHTC/HTC 的差异分析结果为空")
  sc_union_genes <- character(0)
} else {
  sc_union_genes <- unique(deg_sc_all$gene)
}

write.csv(deg_preHTC, file.path(out_dir, "scRNA_DEG_preHTC_OA_vs_HC.csv"), row.names = FALSE)
write.csv(deg_HTC,    file.path(out_dir, "scRNA_DEG_HTC_OA_vs_HC.csv"), row.names = FALSE)
write.csv(data.frame(gene = sc_union_genes),
          file.path(out_dir, "scRNA_DEG_union_preHTC_HTC.csv"),
          row.names = FALSE)

ko_intersect_sc <- intersect(ko_genes, sc_union_genes)

write.csv(
  data.frame(KO_gene = ko_intersect_sc),
  file.path(out_dir, "KO_intersection_with_scRNA_union.csv"),
  row.names = FALSE
)

bulk_dt <- fread(bulk_path, data.table = FALSE)

gene_col <- colnames(bulk_dt)[1]
colnames(bulk_dt)[1] <- "gene"

bulk_dt$gene <- as.character(bulk_dt$gene)
bulk_dt <- bulk_dt[!is.na(bulk_dt$gene) & bulk_dt$gene != "", ]

expr_mat <- as.matrix(bulk_dt[, -1, drop = FALSE])
rownames(expr_mat) <- bulk_dt$gene
mode(expr_mat) <- "numeric"

sample_names <- colnames(expr_mat)

hc_as_cols <- grep("^Control\\.AS", sample_names, value = TRUE)
oa_as_cols <- grep("^(NWB\\.AS|WB\\.AS)", sample_names, value = TRUE)

cat("\nHC.AS 样本数: ", length(hc_as_cols), "\n", sep = "")
cat("OA.AS 样本数: ", length(oa_as_cols), "\n", sep = "")

if (length(hc_as_cols) == 0) stop("未识别到 HC.AS 样本（Control.AS）")
if (length(oa_as_cols) == 0) stop("未识别到 OA.AS 样本（NWB.AS / WB.AS）")

expr_as <- expr_mat[, c(hc_as_cols, oa_as_cols), drop = FALSE]

group_bulk <- factor(
  c(rep("HC_AS", length(hc_as_cols)), rep("OA_AS", length(oa_as_cols))),
  levels = c("HC_AS", "OA_AS")
)

dge <- DGEList(counts = expr_as, group = group_bulk)
keep <- filterByExpr(dge, group = group_bulk)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

design <- model.matrix(~ group_bulk)
colnames(design) <- make.names(colnames(design))

v <- voom(dge, design = design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

coef_name <- grep("group_bulkOA_AS", colnames(coef(fit)), value = TRUE)
if (length(coef_name) != 1) {
  stop("无法唯一识别 OA_AS vs HC_AS 的系数列")
}

deg_bulk <- topTable(
  fit,
  coef = coef_name,
  number = Inf,
  sort.by = "P"
) %>%
  rownames_to_column("gene") %>%
  mutate(
    change = case_when(
      adj.P.Val < bulk_padj_cutoff & logFC >=  bulk_logfc_cutoff ~ "Up_in_OA_AS",
      adj.P.Val < bulk_padj_cutoff & logFC <= -bulk_logfc_cutoff ~ "Down_in_OA_AS",
      TRUE ~ "NS"
    )
  )

deg_bulk_sig <- deg_bulk %>%
  filter(adj.P.Val < bulk_padj_cutoff, abs(logFC) >= bulk_logfc_cutoff)

bulk_deg_genes <- unique(deg_bulk_sig$gene)

write.csv(deg_bulk,     file.path(out_dir, "bulk_DEG_HC_AS_vs_OA_AS_all.csv"), row.names = FALSE)
write.csv(deg_bulk_sig, file.path(out_dir, "bulk_DEG_HC_AS_vs_OA_AS_sig.csv"), row.names = FALSE)
write.csv(data.frame(gene = bulk_deg_genes),
          file.path(out_dir, "bulk_DEG_HC_AS_vs_OA_AS_sig_genes.csv"),
          row.names = FALSE)

ko_intersect_bulk <- intersect(ko_genes, bulk_deg_genes)

write.csv(
  data.frame(KO_gene = ko_intersect_bulk),
  file.path(out_dir, "KO_intersection_with_bulk_AS_DEG.csv"),
  row.names = FALSE
)

summary_df <- data.frame(
  KO_gene = ko_genes,
  in_scRNA_union = ko_genes %in% sc_union_genes,
  in_bulk_AS_DEG = ko_genes %in% bulk_deg_genes,
  in_both = (ko_genes %in% sc_union_genes) & (ko_genes %in% bulk_deg_genes),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(in_both), desc(in_scRNA_union), desc(in_bulk_AS_DEG), KO_gene)

write.csv(summary_df, file.path(out_dir, "KO_gene_intersection_summary.csv"), row.names = FALSE)

ko_sc_detail <- deg_sc_all %>%
  filter(gene %in% ko_genes) %>%
  arrange(cell_type, p_val_adj, desc(abs(logFC_use)))

ko_bulk_detail <- deg_bulk_sig %>%
  filter(gene %in% ko_genes) %>%
  arrange(adj.P.Val, desc(abs(logFC)))

write.csv(ko_sc_detail,   file.path(out_dir, "KO_gene_in_scRNA_DEG_detail.csv"), row.names = FALSE)
write.csv(ko_bulk_detail, file.path(out_dir, "KO_gene_in_bulk_DEG_detail.csv"), row.names = FALSE)

if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
library(VennDiagram)
library(grid)

set_ko   <- unique(ko_genes)
set_sc   <- unique(sc_union_genes)
set_bulk <- unique(bulk_deg_genes)

n_ko <- length(set_ko)
n_sc <- length(set_sc)
n_bulk <- length(set_bulk)

n_ko_sc <- length(intersect(set_ko, set_sc))
n_ko_bulk <- length(intersect(set_ko, set_bulk))
n_sc_bulk <- length(intersect(set_sc, set_bulk))
n_all <- length(Reduce(intersect, list(set_ko, set_sc, set_bulk)))

venn_intersections <- list(
  KO_only = setdiff(set_ko, union(set_sc, set_bulk)),
  sc_only = setdiff(set_sc, union(set_ko, set_bulk)),
  bulk_only = setdiff(set_bulk, union(set_ko, set_sc)),
  KO_sc = setdiff(intersect(set_ko, set_sc), set_bulk),
  KO_bulk = setdiff(intersect(set_ko, set_bulk), set_sc),
  sc_bulk = setdiff(intersect(set_sc, set_bulk), set_ko),
  all_three = Reduce(intersect, list(set_ko, set_sc, set_bulk))
)

write.csv(
  stack(venn_intersections),
  file.path(out_dir, "Venn_intersection_genes.csv"),
  row.names = FALSE
)

venn_plot <- draw.triple.venn(
  area1 = n_ko,
  area2 = n_sc,
  area3 = n_bulk,
  n12 = n_ko_sc,
  n13 = n_ko_bulk,
  n23 = n_sc_bulk,
  n123 = n_all,
  category = c("KO genes", "scRNA preHTC+HTC DEG union", "bulk HC.AS vs OA.AS DEG"),
  fill = c("#4DBBD5", "#00A087", "#E64B35"),
  alpha = c(0.45, 0.45, 0.45),
  cex = 1.4,
  cat.cex = 1.2,
  cat.col = c("#4DBBD5", "#00A087", "#E64B35"),
  cat.pos = c(-20, 20, 180),
  cat.dist = c(0.06, 0.06, 0.06),
  margin = 0.08,
  scaled = FALSE
)

pdf(file.path(out_dir, "Venn_KO_scRNA_bulk.pdf"), width = 10, height = 8)
grid.newpage()
grid.draw(venn_plot)
dev.off()

png(file.path(out_dir, "Venn_KO_scRNA_bulk.png"), width = 3000, height = 2400, res = 300)
grid.newpage()
grid.draw(venn_plot)
dev.off()

cat("\n韦恩图已输出：\n")
cat(file.path(out_dir, "Venn_KO_scRNA_bulk.pdf"), "\n")
cat(file.path(out_dir, "Venn_KO_scRNA_bulk.png"), "\n")
cat(file.path(out_dir, "Venn_intersection_genes.csv"), "\n")

cat("\n================ 结果汇总 ================\n")
cat("KO gene 总数: ", length(ko_genes), "\n", sep = "")
cat("scRNA preHTC+HTC DEG并集基因数: ", length(sc_union_genes), "\n", sep = "")
cat("bulk HC.AS vs OA.AS DEG基因数: ", length(bulk_deg_genes), "\n", sep = "")

cat("\n与 scRNA DEG并集 的交集基因：\n")
print(ko_intersect_sc)

cat("\n与 bulk HC.AS vs OA.AS DEG 的交集基因：\n")
print(ko_intersect_bulk)

cat("\n同时出现在 scRNA DEG并集 和 bulk AS DEG 中的 KO genes：\n")
print(summary_df$KO_gene[summary_df$in_both])

cat("\n输出目录：\n", out_dir, "\n", sep = "")
