rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
options(stringsAsFactors = FALSE)

WORKDIR <- project_path("results_root", "core_genes")
dir.create(WORKDIR, showWarnings = FALSE, recursive = TRUE)
setwd(WORKDIR)

COUNT_FILE <- project_path("data_root", "bulk_reference/OA_COUNTS_raw.csv")
META_FILE  <- project_path("data_root", "bulk_reference/sample_info.csv")

Seed_GlyLE <- c(
  "MDH1","DCN","PKM","LDHA","SLC35A3","PLOD2","CD44","ENO1","PPIA","AK3",
  "TGFBI","AGL","NASP","LHX9","COPB2","PAM","VCAN","ELF3","SDC2","ME2","PGK1","HDLBP","COG2"
)

FDR_AXIS  <- 0.05
RHO_MIN   <- 0.30
FDR_DE    <- 0.05
LOGFC_MIN <- 0.25

LAYER_AS  <- "AS"
GROUP_OA  <- "OA"
GROUP_HC  <- "HC"

pkgs <- c("data.table", "edgeR", "limma", "ggplot2")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
library(data.table)
library(edgeR)
library(limma)
library(ggplot2)

normalize_group <- function(x) {
  x2 <- toupper(trimws(as.character(x)))
  x2 <- gsub("\\s+", "", x2)

  x2[x2 %in% c("CONTROL","CTRL","CON","HEALTHY","NORMAL")] <- "HC"
  x2[x2 %in% c("OSTEOARTHRITIS")] <- "OA"
  return(x2)
}

normalize_layer <- function(x) {
  x2 <- toupper(trimws(as.character(x)))
  x2 <- gsub("\\s+", "", x2)
  return(x2)
}

rank_desc <- function(x) rank(x, ties.method = "average")

message("Reading counts: ", COUNT_FILE)
cnt_dt <- fread(COUNT_FILE)
stopifnot(ncol(cnt_dt) >= 2)

gene_col <- names(cnt_dt)[1]
genes <- cnt_dt[[gene_col]]
count_mat <- as.matrix(cnt_dt[, -1, with = FALSE])
rownames(count_mat) <- genes
storage.mode(count_mat) <- "numeric"
if (anyNA(count_mat)) stop("Counts contain NA. Check OA_COUNTS_raw.csv.")

samples <- colnames(count_mat)
message("Counts loaded | Samples: ", length(samples), " | Genes: ", nrow(count_mat))

message("Reading meta: ", META_FILE)
meta <- fread(META_FILE)
meta <- as.data.frame(meta)

need_cols <- c("sample_id", "condition", "region", "group")
miss <- setdiff(need_cols, colnames(meta))
if (length(miss) > 0) stop("sample_info.csv 缺少列：", paste(miss, collapse = ", "))

meta$sample_id <- as.character(meta$sample_id)
meta$condition <- normalize_group(meta$condition)
meta$region    <- normalize_layer(meta$region)
meta$group     <- as.character(meta$group)

meta_in <- meta[meta$sample_id %in% samples, , drop = FALSE]
missing_meta <- setdiff(samples, meta_in$sample_id)
if (length(missing_meta) > 0) {
  stop("这些样本在 counts 里，但在 sample_info.csv 找不到：\n", paste(missing_meta, collapse = ", "))
}
meta_in <- meta_in[match(samples, meta_in$sample_id), , drop = FALSE]
stopifnot(all(meta_in$sample_id == samples))

fwrite(meta_in, file.path(WORKDIR, "meta_used.tsv"), sep = "\t")
message("Wrote meta_used.tsv")

Seed_GlyLE_in <- intersect(Seed_GlyLE, rownames(count_mat))
if (length(Seed_GlyLE_in) == 0) {
  stop("None of Seed_GlyLE genes found in counts rownames. Check gene symbols/case.")
}
write.table(Seed_GlyLE_in, "Seed_GlyLE_in_counts.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

idx_as <- which(meta_in$region == LAYER_AS)
if (length(idx_as) < 4) stop("AS 样本数 < 4。请检查 sample_info.csv 的 region 列是否正确标注 AS。")

meta_as <- meta_in[idx_as, , drop = FALSE]
cnt_as  <- count_mat[, idx_as, drop = FALSE]

ugrp <- unique(meta_as$condition)
if (!all(c(GROUP_OA, GROUP_HC) %in% ugrp)) {
  stop("AS 样本中没有同时包含 OA 与 HC。当前 AS condition 为：",
       paste(ugrp, collapse = ", "),
       "\n请检查 sample_info.csv 的 condition 列。")
}

dge <- DGEList(counts = cnt_as)
dge <- calcNormFactors(dge)
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

message("Building axis proxy: PCA(PC1) on AS logCPM ...")
pc <- prcomp(t(logcpm), center = TRUE, scale. = FALSE)
axis_proxy <- pc$x[, 1]

if (mean(axis_proxy[meta_as$condition == GROUP_OA]) < mean(axis_proxy[meta_as$condition == GROUP_HC])) {
  axis_proxy <- -axis_proxy
}
meta_as$axis <- axis_proxy

meta_in$axis <- NA_real_
meta_in$axis[idx_as] <- meta_as$axis
fwrite(meta_in, file.path(WORKDIR, "meta_used_with_axis_proxy.tsv"), sep = "\t")

grp <- factor(meta_as$condition, levels = c(GROUP_HC, GROUP_OA))
design <- model.matrix(~ 0 + grp)
colnames(design) <- levels(grp)

v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
contr <- makeContrasts(OA_vs_HC = OA - HC, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contr))

de_tab <- topTable(fit2, coef = "OA_vs_HC", number = Inf, sort.by = "none")
de_tab$gene <- rownames(de_tab)
de_tab <- de_tab[, c("gene", "logFC", "P.Value", "adj.P.Val", "AveExpr", "t", "B")]
colnames(de_tab) <- c("gene", "log2FC", "pval", "FDR", "AveExpr", "t", "B")
fwrite(de_tab, "DE_AS_OA_vs_HC_all_genes.tsv", sep = "\t")

DE_AS <- de_tab$gene[de_tab$FDR < FDR_DE & abs(de_tab$log2FC) >= LOGFC_MIN]

axis_vec <- meta_as$axis
message("Computing Axis association (Spearman) within AS ...")

rho <- rep(NA_real_, nrow(logcpm))
pval <- rep(NA_real_, nrow(logcpm))

for (i in seq_len(nrow(logcpm))) {
  ct <- suppressWarnings(cor.test(as.numeric(logcpm[i, ]), axis_vec,
                                  method = "spearman", exact = FALSE))
  rho[i]  <- unname(ct$estimate)
  pval[i] <- ct$p.value
}
padj <- p.adjust(pval, method = "BH")

axis_tab <- data.frame(
  gene = rownames(logcpm),
  rho  = rho,
  pval = pval,
  FDR  = padj,
  stringsAsFactors = FALSE
)
fwrite(axis_tab, "AxisAssoc_AS_spearman_all_genes.tsv", sep = "\t")

AxisAssoc <- axis_tab$gene[axis_tab$FDR < FDR_AXIS & abs(axis_tab$rho) >= RHO_MIN]

CoreCandidates_Gly <- Reduce(intersect, list(Seed_GlyLE_in, AxisAssoc, DE_AS))

core_merge <- merge(axis_tab, de_tab, by = "gene", all = FALSE)
core_merge$in_Seed_GlyLE <- core_merge$gene %in% Seed_GlyLE_in
core_merge$in_AxisAssoc  <- core_merge$gene %in% AxisAssoc
core_merge$in_DE_AS      <- core_merge$gene %in% DE_AS
core_merge$in_CoreCandidates_Gly <- core_merge$gene %in% CoreCandidates_Gly

core_tab <- core_merge[core_merge$in_CoreCandidates_Gly, , drop = FALSE]
core_tab$direction_ok_gly_up <- (core_tab$rho > 0 & core_tab$log2FC > 0)

core_tab$Score <- rank_desc(abs(core_tab$rho)) + rank_desc(abs(core_tab$log2FC))

core_gly_up <- core_tab[core_tab$direction_ok_gly_up, , drop = FALSE]
core_gly_up <- core_gly_up[order(core_gly_up$Score, decreasing = TRUE), , drop = FALSE]

fwrite(core_tab[order(core_tab$Score, decreasing = TRUE), ],
       "CoreCandidates_Gly_all_core_with_scores.tsv", sep = "\t")
fwrite(core_gly_up, "Core_Gly_Up_ranked.tsv", sep = "\t")
write.table(core_gly_up$gene, "Core_Gly_Up_ranked_genes.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

sizes <- data.frame(
  set = c("Seed_GlyLE_in", "AxisAssoc", "DE_AS", "CoreCandidates_Gly", "Core_Gly_Up"),
  n   = c(length(Seed_GlyLE_in), length(AxisAssoc), length(DE_AS),
          length(CoreCandidates_Gly), nrow(core_gly_up))
)
fwrite(sizes, "Set_sizes.tsv", sep = "\t")

if (nrow(core_tab) > 0) {
  p <- ggplot(core_tab, aes(x = rho, y = log2FC)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_point(aes(color = direction_ok_gly_up), alpha = 0.85) +
    theme_bw() +
    labs(title = "CoreCandidates_Gly (AS): rho vs log2FC",
         x = "Spearman rho (logCPM vs axis proxy)",
         y = "log2FC (OA vs HC within AS)")
  ggsave("CoreCandidates_rho_vs_log2FC.png", p, width = 6.5, height = 5.2, dpi = 220)
}

message("DONE. Workdir: ", WORKDIR)
message("Key outputs:")
message(" - meta_used.tsv")
message(" - meta_used_with_axis_proxy.tsv")
message(" - DE_AS_OA_vs_HC_all_genes.tsv")
message(" - AxisAssoc_AS_spearman_all_genes.tsv")
message(" - Core_Gly_Up_ranked.tsv")
message(" - Core_Gly_Up_ranked_genes.txt")
