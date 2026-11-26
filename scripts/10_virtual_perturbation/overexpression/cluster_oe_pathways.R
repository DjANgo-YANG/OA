rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
options(stringsAsFactors = FALSE)

cran_pkgs <- c(
  "data.table", "dplyr", "stringr", "purrr", "tidyr",
  "ggplot2", "scales", "RColorBrewer", "pheatmap", "tibble", "cowplot"
)

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

if (!requireNamespace("uwot", quietly = TRUE)) {
  install.packages("uwot")
}

library(data.table)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(pheatmap)
library(tibble)
library(cowplot)
library(uwot)

base_dir <- project_path("results_root", "virtual_oe")
gsea_base_dir <- file.path(base_dir, "HALLMARK_GSEA_EnrichmentMap", "per_OE_results")
out_dir <- file.path(base_dir, "OE_Clustering_Module_Overview_clusterEllipse")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pathway_fdr_cutoff <- 0.05
top_n_pathways_global <- 20
top_n_dotplot <- 15
min_kos_per_pathway <- 2
use_absNES <- FALSE
cluster_k <- 4

pdf_width <- 14
pdf_height <- 10

extract_oe_gene <- function(folder_name) {
  folder_name %>%
    str_remove("^\\d+_") %>%
    str_remove("过表达$") %>%
    str_split("_") %>%
    map_chr(~ .x[1])
}

read_gsea_sig <- function(gsea_sig_file, oe_gene,
                          pathway_fdr_cutoff = 0.05) {
  dt <- fread(gsea_sig_file)

  if (!"Description" %in% colnames(dt) && "pathway" %in% colnames(dt)) {
    dt$Description <- dt$pathway
  }
  if (!"FDR" %in% colnames(dt) && "padj" %in% colnames(dt)) {
    dt$FDR <- dt$padj
  }

  need_cols <- c("Description", "FDR")
  miss_cols <- setdiff(need_cols, colnames(dt))
  if (length(miss_cols) > 0) {
    stop("GSEA sig 文件缺少必要列: ",
         paste(miss_cols, collapse = ", "),
         " | 文件: ", gsea_sig_file)
  }

  dt2 <- dt %>%
    mutate(
      OE_gene = oe_gene,
      pathway = as.character(Description),
      FDR = as.numeric(FDR),
      NES = if ("NES" %in% colnames(.)) as.numeric(NES) else NA_real_
    ) %>%
    filter(!is.na(FDR), FDR <= pathway_fdr_cutoff)

  dt2
}

pick_module_name <- function(pathway_vec) {
  p <- toupper(pathway_vec)

  if (any(str_detect(p, "GLYCOLYSIS|OXIDATIVE_PHOSPHORYLATION|FATTY_ACID|CHOLESTEROL|BILE_ACID|XENOBIOTIC|MTORC1|MYC_TARGETS|PEROXISOME|PROTEIN_SECRETION|HEME_METABOLISM"))) {
    return("Metabolic")
  }
  if (any(str_detect(p, "EPITHELIAL_MESENCHYMAL_TRANSITION|ANGIOGENESIS|APICAL_JUNCTION|APICAL_SURFACE|MYOGENESIS|TGF_BETA_SIGNALING"))) {
    return("ECM_Remodeling")
  }
  if (any(str_detect(p, "HYPOXIA|REACTIVE_OXYGEN_SPECIES|UNFOLDED_PROTEIN_RESPONSE|P53_PATHWAY|DNA_REPAIR|UV_RESPONSE|APOPTOSIS"))) {
    return("Stress_Response")
  }
  if (any(str_detect(p, "TNFA_SIGNALING_VIA_NFKB|INFLAMMATORY_RESPONSE|INTERFERON|IL2_STAT5|IL6_JAK_STAT3|COMPLEMENT|ALLOGRAFT_REJECTION"))) {
    return("Inflammatory")
  }
  if (any(str_detect(p, "E2F_TARGETS|G2M_CHECKPOINT|MITOTIC_SPINDLE"))) {
    return("Cell_Cycle")
  }
  return("Mixed")
}

oe_dirs <- list.dirs(gsea_base_dir, recursive = FALSE, full.names = TRUE)

if (length(oe_dirs) == 0) {
  stop("未找到 per_OE_results 下的 OE 文件夹，请检查路径。")
}

message("找到 OE 结果文件夹数: ", length(oe_dirs))

gsea_all_list <- list()
skip_log <- list()

for (d in oe_dirs) {
  folder_name <- basename(d)
  oe_gene <- extract_oe_gene(folder_name)

  gsea_sig_file <- file.path(
    d,
    paste0(folder_name, "_HALLMARK_GSEA_sig.csv")
  )

  if (!file.exists(gsea_sig_file)) {
    skip_log[[folder_name]] <- "missing HALLMARK_GSEA_sig.csv"
    next
  }

  message("读取: ", folder_name)

  gsea_one <- tryCatch(
    read_gsea_sig(
      gsea_sig_file = gsea_sig_file,
      oe_gene = oe_gene,
      pathway_fdr_cutoff = pathway_fdr_cutoff
    ),
    error = function(e) {
      skip_log[[folder_name]] <<- paste("read error:", e$message)
      return(NULL)
    }
  )

  if (is.null(gsea_one) || nrow(gsea_one) == 0) {
    skip_log[[folder_name]] <- "empty significant pathways"
    next
  }

  gsea_all_list[[folder_name]] <- gsea_one
}

gsea_all <- bind_rows(gsea_all_list)

if (nrow(gsea_all) == 0) {
  stop("没有读到有效的 OE GSEA 结果。")
}

if (length(skip_log) > 0) {
  skip_df <- data.frame(
    folder = names(skip_log),
    reason = unlist(skip_log),
    row.names = NULL
  )
  fwrite(skip_df, file.path(out_dir, "OE_skipped_reason.csv"))
}

pathway_summary <- gsea_all %>%
  group_by(pathway) %>%
  summarise(
    n_OEs = n_distinct(OE_gene),
    mean_absNES = mean(abs(NES), na.rm = TRUE),
    minFDR = min(FDR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_OEs >= min_kos_per_pathway) %>%
  arrange(desc(n_OEs), desc(mean_absNES), minFDR)

selected_pathways <- pathway_summary %>%
  slice_head(n = top_n_pathways_global) %>%
  pull(pathway)

if (length(selected_pathways) == 0) {
  stop("没有 pathway 满足 min_kos_per_pathway 条件，请降低 min_kos_per_pathway 或放宽筛选。")
}

heat_df <- gsea_all %>%
  filter(pathway %in% selected_pathways) %>%
  transmute(
    OE_gene = OE_gene,
    pathway = pathway,
    NES = as.numeric(NES),
    absNES = abs(as.numeric(NES)),
    padj = as.numeric(FDR)
  )

oe_pathway_mat <- heat_df %>%
  select(OE_gene, pathway, absNES) %>%
  pivot_wider(names_from = pathway, values_from = absNES, values_fill = 0)

mat <- as.data.frame(oe_pathway_mat)
rownames(mat) <- mat$OE_gene
mat$OE_gene <- NULL
mat <- as.matrix(mat)

col_var <- apply(mat, 2, var, na.rm = TRUE)
col_var[is.na(col_var)] <- 0
keep_cols <- col_var > 0

removed_pathways <- names(col_var)[!keep_cols]
write.csv(
  data.frame(pathway = removed_pathways),
  file.path(out_dir, "removed_constant_pathways.csv"),
  row.names = FALSE
)

mat <- mat[, keep_cols, drop = FALSE]

cat("原始 pathway 数: ", length(col_var), "\n", sep = "")
cat("保留 pathway 数: ", ncol(mat), "\n", sep = "")
cat("删除常数列数: ", sum(!keep_cols), "\n", sep = "")

if (ncol(mat) < 2) {
  stop("过滤常数列后剩余 pathway 数量不足，无法进行热图绘制。")
}

fwrite(
  as.data.table(cbind(OE_gene = rownames(mat), as.data.frame(mat))),
  file.path(out_dir, "OE_pathway_signature_matrix.csv")
)

hc_rows <- hclust(dist(mat), method = "ward.D2")
row_clusters <- cutree(hc_rows, k = min(cluster_k, nrow(mat)))

module_anno <- data.frame(
  OE_gene = rownames(mat),
  cluster = paste0("Cluster_", row_clusters),
  stringsAsFactors = FALSE
)

fwrite(module_anno, file.path(out_dir, "OE_module_annotation.csv"))

oe_order <- hc_rows$labels[hc_rows$order]

hc_cols <- hclust(dist(t(mat)), method = "ward.D2")
pathway_order <- hc_cols$labels[hc_cols$order]

heat_plot_df <- expand.grid(
  OE_gene = oe_order,
  pathway = pathway_order,
  stringsAsFactors = FALSE
) %>%
  left_join(
    heat_df %>% select(OE_gene, pathway, absNES, padj),
    by = c("OE_gene", "pathway")
  ) %>%
  mutate(
    absNES = ifelse(is.na(absNES), 0, absNES),
    sig_label = case_when(
      is.na(padj) ~ "",
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      padj < 0.10  ~ "·",
      TRUE ~ ""
    )
  ) %>%
  left_join(module_anno, by = "OE_gene")

score_rng <- range(heat_plot_df$absNES, na.rm = TRUE)

heat_plot_df <- heat_plot_df %>%
  mutate(
    plot_score = if (diff(score_rng) == 0) 0 else
      scales::rescale(absNES, to = c(0, 2), from = score_rng),
    OE_gene = factor(OE_gene, levels = oe_order),
    pathway = factor(pathway, levels = rev(pathway_order))
  )

cluster_levels <- unique(module_anno$cluster)
cluster_cols <- setNames(
  RColorBrewer::brewer.pal(max(3, min(8, length(cluster_levels))), "Set2")[seq_along(cluster_levels)],
  cluster_levels
)

anno_top_df <- data.frame(
  OE_gene = oe_order,
  stringsAsFactors = FALSE
) %>%
  left_join(
    module_anno %>% select(OE_gene, cluster),
    by = "OE_gene"
  ) %>%
  mutate(
    OE_gene = factor(OE_gene, levels = oe_order),
    y1 = "cluster"
  )

p_heat <- ggplot(heat_plot_df, aes(x = OE_gene, y = pathway, fill = plot_score)) +
  geom_tile(color = "grey85", linewidth = 0.35) +
  geom_text(aes(label = sig_label), size = 3.0, color = "black") +
  scale_fill_gradient(
    low = "#ffffff",
    high = "#0c5096",
    limits = c(0, 2),
    oob = scales::squish,
    name = "abs(NES)"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(
    title = "OE clustering based on pathway signatures",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )

p_anno_cluster <- ggplot(anno_top_df, aes(x = OE_gene, y = y1, fill = cluster)) +
  geom_tile(color = "grey85", linewidth = 0.3) +
  scale_fill_manual(values = cluster_cols, name = "cluster", drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    plot.margin = margin(0, 5, 0, 5)
  )

legend_cluster <- cowplot::get_legend(p_anno_cluster)
p_anno_cluster_noleg <- p_anno_cluster + theme(legend.position = "none")

top_annos <- cowplot::plot_grid(
  p_anno_cluster_noleg,
  ncol = 1,
  rel_heights = c(1),
  align = "v",
  axis = "lr"
)

main_left <- cowplot::plot_grid(
  top_annos,
  p_heat,
  ncol = 1,
  rel_heights = c(0.08, 1),
  align = "v",
  axis = "lr"
)

p_heat_final <- cowplot::plot_grid(
  main_left,
  legend_cluster,
  nrow = 1,
  rel_widths = c(1, 0.12)
)

ggsave(
  file.path(out_dir, "OE_pathway_heatmap_referenceOldCode.pdf"),
  p_heat_final,
  width = 18,
  height = 10
)

ggsave(
  file.path(out_dir, "OE_pathway_heatmap_referenceOldCode.png"),
  p_heat_final,
  width = 18,
  height = 10,
  dpi = 300
)

heat_rot_df <- heat_plot_df %>%
  mutate(
    OE_gene = factor(OE_gene, levels = rev(oe_order)),
    pathway = factor(pathway, levels = rev(pathway_order))
  )

p_heat_rot <- ggplot(heat_rot_df, aes(x = pathway, y = OE_gene, fill = plot_score)) +
  geom_tile(color = "grey85", linewidth = 0.35) +
  geom_text(aes(label = sig_label), size = 3.0, color = "black") +
  scale_fill_gradient(
    low = "#ffffff",
    high = "#0c5096",
    limits = c(0, 2),
    oob = scales::squish,
    name = "abs(NES)"
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(
    title = "OE clustering based on pathway signatures (rotated)",
    x = "Pathway",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )

anno_side_df <- anno_top_df %>%
  mutate(OE_gene = factor(OE_gene, levels = rev(oe_order)))

p_anno_side <- ggplot(anno_side_df, aes(x = 1, y = OE_gene, fill = cluster)) +
  geom_tile(color = "grey85", linewidth = 0.3) +
  scale_fill_manual(values = cluster_cols, name = "cluster", drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  theme_void(base_size = 12) +
  theme(legend.position = "none")

p_rot_final <- cowplot::plot_grid(
  p_heat_rot,
  p_anno_side,
  nrow = 1,
  rel_widths = c(1, 0.025),
  align = "h",
  axis = "tb"
)

ggsave(
  file.path(out_dir, "OE_pathway_heatmap_rotated_sig_blue_fixed.pdf"),
  p_rot_final,
  width = 18,
  height = 10
)

ggsave(
  file.path(out_dir, "OE_pathway_heatmap_rotated_sig_blue_fixed.png"),
  p_rot_final,
  width = 18,
  height = 10,
  dpi = 300
)

pca_res <- prcomp(mat, center = TRUE, scale. = TRUE)

pca_df <- data.frame(
  OE_gene = rownames(mat),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  stringsAsFactors = FALSE
) %>%
  left_join(module_anno, by = "OE_gene")

ellipse_groups_pca <- pca_df %>%
  count(cluster) %>%
  filter(n >= 3) %>%
  pull(cluster)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = cluster, label = OE_gene)) +
  stat_ellipse(
    data = subset(pca_df, cluster %in% ellipse_groups_pca),
    aes(group = cluster),
    linewidth = 0.9,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  geom_point(size = 4) +
  geom_text(nudge_y = 0.15, size = 3.8) +
  scale_color_manual(values = cluster_cols, drop = FALSE) +
  theme_bw(base_size = 12) +
  labs(
    title = "PCA of OE pathway signatures",
    x = paste0("PC1 (", round(100 * summary(pca_res)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca_res)$importance[2, 2], 1), "%)"),
    color = "Cluster"
  ) +
  theme(
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "OE_PCA_clusterEllipse.pdf"), p_pca, width = 10, height = 8)
ggsave(file.path(out_dir, "OE_PCA_clusterEllipse.png"), p_pca, width = 10, height = 8, dpi = 300)

set.seed(123)

n_neighbors_use <- min(8, nrow(mat) - 1)
if (n_neighbors_use < 2) {
  warning("OE 数量过少，无法稳定运行 UMAP，跳过 UMAP 输出。")
} else {
  umap_res <- uwot::umap(
    mat,
    n_neighbors = n_neighbors_use,
    min_dist = 0.35,
    metric = "euclidean",
    scale = TRUE,
    verbose = TRUE
  )

  umap_df <- data.frame(
    OE_gene = rownames(mat),
    UMAP1 = umap_res[, 1],
    UMAP2 = umap_res[, 2],
    stringsAsFactors = FALSE
  ) %>%
    left_join(module_anno, by = "OE_gene")

  ellipse_groups_umap <- umap_df %>%
    count(cluster) %>%
    filter(n >= 3) %>%
    pull(cluster)

  p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, color = cluster, label = OE_gene)) +
    stat_ellipse(
      data = subset(umap_df, cluster %in% ellipse_groups_umap),
      aes(group = cluster),
      linewidth = 0.9,
      alpha = 0.9,
      show.legend = FALSE
    ) +
    geom_point(size = 4) +
    geom_text(nudge_y = 0.12, size = 3.8) +
    scale_color_manual(values = cluster_cols, drop = FALSE) +
    theme_bw(base_size = 12) +
    labs(
      title = "UMAP of OE pathway signatures",
      x = "UMAP1",
      y = "UMAP2",
      color = "Cluster"
    ) +
    theme(
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(out_dir, "OE_UMAP_clusterEllipse.pdf"), p_umap, width = 10, height = 8)
  ggsave(file.path(out_dir, "OE_UMAP_clusterEllipse.png"), p_umap, width = 10, height = 8, dpi = 300)
}

module_label_df <- gsea_all %>%
  group_by(OE_gene) %>%
  summarise(
    module_label = pick_module_name(unique(pathway)),
    .groups = "drop"
  )

pathway_dot_summary <- gsea_all %>%
  left_join(module_label_df, by = "OE_gene") %>%
  group_by(module_label, pathway) %>%
  summarise(
    meanNES = mean(NES, na.rm = TRUE),
    nOE = n_distinct(OE_gene),
    .groups = "drop"
  ) %>%
  arrange(module_label, desc(abs(meanNES)), desc(nOE))

pathway_rank_df <- pathway_dot_summary %>%
  group_by(pathway) %>%
  summarise(
    total_mean_absNES = mean(abs(meanNES), na.rm = TRUE),
    total_nOE = sum(nOE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_nOE), desc(total_mean_absNES))

if (nrow(pathway_rank_df) == 0) {
  warning("pathway_rank_df 为空，跳过模块 dotplot 输出。")
} else {
  n_top_dot <- min(top_n_dotplot, nrow(pathway_rank_df))

  top_pathways_dot <- pathway_rank_df %>%
    slice_head(n = n_top_dot) %>%
    pull(pathway)

  dot_df <- pathway_dot_summary %>%
    filter(pathway %in% top_pathways_dot)

  pathway_levels_dot <- dot_df %>%
    group_by(pathway) %>%
    summarise(
      total_nOE = sum(nOE),
      total_mean_absNES = mean(abs(meanNES), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(total_nOE, total_mean_absNES) %>%
    pull(pathway)

  module_levels_dot <- c("Metabolic", "ECM_Remodeling", "Stress_Response", "Inflammatory", "Cell_Cycle", "Mixed")
  module_levels_dot <- module_levels_dot[module_levels_dot %in% unique(dot_df$module_label)]

  dot_df <- dot_df %>%
    mutate(
      pathway = factor(pathway, levels = pathway_levels_dot),
      module_label = factor(module_label, levels = module_levels_dot)
    )

  p_dot <- ggplot(dot_df, aes(x = module_label, y = pathway)) +
    geom_point(aes(size = nOE, color = meanNES)) +
    scale_size_continuous(range = c(3, 10), name = "OE count") +
    scale_color_gradient2(
      low = "#313695",
      mid = "#f7f7f7",
      high = "#a50026",
      midpoint = 0,
      name = "Mean NES"
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = "Functional module dotplot of OE groups",
      x = "OE module",
      y = NULL
    ) +
    theme(
      axis.text.x = element_text(angle = 20, hjust = 1),
      panel.grid.major = element_line(color = "grey92"),
      panel.grid.minor = element_blank()
    )

  ggsave(file.path(out_dir, "OE_module_dotplot.pdf"), p_dot, width = 10, height = 8)
  ggsave(file.path(out_dir, "OE_module_dotplot.png"), p_dot, width = 10, height = 8, dpi = 300)
}

module_summary <- module_label_df %>%
  group_by(module_label) %>%
  summarise(
    OE_genes = paste(OE_gene, collapse = "; "),
    n = n(),
    .groups = "drop"
  )

fwrite(module_summary, file.path(out_dir, "OE_module_summary.csv"))

cat(
  "输出完成：\n",
  file.path(out_dir, "OE_pathway_signature_matrix.csv"), "\n",
  file.path(out_dir, "removed_constant_pathways.csv"), "\n",
  file.path(out_dir, "OE_module_annotation.csv"), "\n",
  file.path(out_dir, "OE_module_summary.csv"), "\n",
  file.path(out_dir, "OE_pathway_heatmap_referenceOldCode.pdf"), "\n",
  file.path(out_dir, "OE_pathway_heatmap_rotated_sig_blue_fixed.pdf"), "\n",
  file.path(out_dir, "OE_PCA_clusterEllipse.pdf"), "\n",
  file.path(out_dir, "OE_module_dotplot.pdf"), "\n",
  if (file.exists(file.path(out_dir, "OE_UMAP_clusterEllipse.pdf"))) file.path(out_dir, "OE_UMAP_clusterEllipse.pdf") else "UMAP skipped", "\n"
)
