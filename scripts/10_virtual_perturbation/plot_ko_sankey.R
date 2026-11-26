rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
options(stringsAsFactors = FALSE)

cran_pkgs <- c(
  "data.table", "dplyr", "stringr", "purrr", "tidyr",
  "ggplot2", "ggalluvial", "cowplot", "scales"
)

for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(data.table)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)
library(ggalluvial)
library(cowplot)
library(scales)

base_dir <- Sys.getenv("OAGLY_KO_ROOT", unset = project_path("results_root", "virtual_ko"))
gsea_base_dir <- file.path(base_dir, "HALLMARK_GSEA_EnrichmentMap", "per_KO_results")
out_dir <- file.path(base_dir, "KO_Gene_Pathway_Sankey_FINAL_v9")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

top_n_genes_per_ko <- 15
only_significant_genes <- TRUE
exclude_self_ko_gene <- TRUE
min_abs_z <- 0

top_n_pathways_per_ko <- 10
pathway_fdr_cutoff <- 0.05

use_shared_gene_filter <- TRUE
min_shared_kos <- 1

max_gene_labels_total <- 45

show_top_n_pathways <- 10

pdf_width <- 24
pdf_height <- 15

ko_label_size <- 4
gene_label_size <- 4
pathway_label_size <- 4

extract_ko_gene <- function(folder_name) {
  folder_name %>%
    str_remove("^\\d+_") %>%
    str_remove("敲除$")
}

split_genes <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  x %>%
    str_split(",") %>%
    unlist() %>%
    str_trim() %>%
    unique()
}

read_diffreg_topgenes <- function(diff_file, ko_gene,
                                  top_n_genes_per_ko = 12,
                                  only_significant_genes = TRUE,
                                  exclude_self_ko_gene = TRUE,
                                  min_abs_z = 0) {
  dt <- fread(diff_file)

  need_cols <- c("gene", "Z")
  miss_cols <- setdiff(need_cols, colnames(dt))
  if (length(miss_cols) > 0) {
    stop("diffRegulation_full.txt 缺少必要列: ",
         paste(miss_cols, collapse = ", "),
         " | 文件: ", diff_file)
  }

  dt2 <- dt %>%
    mutate(
      gene = as.character(gene),
      Z = as.numeric(Z),
      absZ = abs(Z),
      Significance = if ("Significance" %in% colnames(.)) as.character(Significance) else NA_character_
    ) %>%
    filter(!is.na(gene), gene != "", !is.na(Z))

  if (only_significant_genes && "Significance" %in% colnames(dt2)) {
    dt2 <- dt2 %>% filter(Significance == "Significant")
  }

  if (exclude_self_ko_gene) {
    dt2 <- dt2 %>% filter(gene != ko_gene)
  }

  dt2 %>%
    filter(absZ >= min_abs_z) %>%
    arrange(desc(absZ), desc(Z)) %>%
    slice_head(n = top_n_genes_per_ko) %>%
    mutate(KO_gene = ko_gene)
}

read_gsea_sig <- function(gsea_sig_file, ko_gene,
                          pathway_fdr_cutoff = 0.05,
                          top_n_pathways_per_ko = 8) {
  dt <- fread(gsea_sig_file)

  if (!"Description" %in% colnames(dt) && "pathway" %in% colnames(dt)) {
    dt$Description <- dt$pathway
  }
  if (!"GeneSetID" %in% colnames(dt) && "Description" %in% colnames(dt)) {
    dt$GeneSetID <- dt$Description
  }
  if (!"FDR" %in% colnames(dt) && "padj" %in% colnames(dt)) {
    dt$FDR <- dt$padj
  }
  if (!"leadingEdge_genes" %in% colnames(dt) && "leadingEdge" %in% colnames(dt)) {
    dt$leadingEdge_genes <- dt$leadingEdge
  }

  need_cols <- c("Description", "FDR", "leadingEdge_genes")
  miss_cols <- setdiff(need_cols, colnames(dt))
  if (length(miss_cols) > 0) {
    stop("GSEA sig 文件缺少必要列: ",
         paste(miss_cols, collapse = ", "),
         " | 文件: ", gsea_sig_file)
  }

  dt %>%
    mutate(
      KO_gene = ko_gene,
      FDR = as.numeric(FDR),
      NES = if ("NES" %in% colnames(.)) as.numeric(NES) else NA_real_
    ) %>%
    filter(!is.na(FDR), FDR <= pathway_fdr_cutoff) %>%
    arrange(FDR, desc(abs(NES))) %>%
    slice_head(n = top_n_pathways_per_ko)
}

make_grad_desc <- function(base_col, n, light_col) {
  rev(colorRampPalette(c(light_col, base_col))(n))
}

ko_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
ko_dirs <- ko_dirs[str_detect(basename(ko_dirs), "敲除")]

if (length(ko_dirs) == 0) {
  stop("未找到 KO 文件夹，请检查 base_dir。")
}

message("找到 KO 文件夹数: ", length(ko_dirs))

triplets_list <- list()
skip_log <- list()

for (d in ko_dirs) {
  folder_name <- basename(d)
  ko_gene <- extract_ko_gene(folder_name)

  diff_file <- file.path(d, "diffRegulation_full.txt")
  gsea_sig_file <- file.path(
    gsea_base_dir,
    folder_name,
    paste0(folder_name, "_HALLMARK_GSEA_sig.csv")
  )

  if (!file.exists(diff_file)) {
    skip_log[[folder_name]] <- "missing diffRegulation_full.txt"
    next
  }
  if (!file.exists(gsea_sig_file)) {
    skip_log[[folder_name]] <- "missing GSEA sig file"
    next
  }

  message("处理: ", folder_name)

  top_genes <- read_diffreg_topgenes(
    diff_file = diff_file,
    ko_gene = ko_gene,
    top_n_genes_per_ko = top_n_genes_per_ko,
    only_significant_genes = only_significant_genes,
    exclude_self_ko_gene = exclude_self_ko_gene,
    min_abs_z = min_abs_z
  )

  if (nrow(top_genes) == 0) {
    skip_log[[folder_name]] <- "top_genes empty after filtering"
    next
  }

  gsea_sig <- read_gsea_sig(
    gsea_sig_file = gsea_sig_file,
    ko_gene = ko_gene,
    pathway_fdr_cutoff = pathway_fdr_cutoff,
    top_n_pathways_per_ko = top_n_pathways_per_ko
  )

  if (nrow(gsea_sig) == 0) {
    skip_log[[folder_name]] <- "gsea_sig empty after filtering"
    next
  }

  gene_set <- unique(top_genes$gene)

  tri_one <- map_dfr(seq_len(nrow(gsea_sig)), function(i) {
    pw <- gsea_sig[i, ]
    le_genes <- split_genes(pw$leadingEdge_genes)
    overlap_genes <- intersect(gene_set, le_genes)

    if (length(overlap_genes) == 0) return(NULL)

    gene_tbl <- top_genes %>%
      filter(gene %in% overlap_genes) %>%
      select(gene, Z, absZ)

    tibble(
      KO_gene = ko_gene,
      affected_gene = gene_tbl$gene,
      pathway = pw$Description,
      gene_absZ = gene_tbl$absZ,
      gene_Z = gene_tbl$Z,
      pathway_FDR = pw$FDR,
      pathway_NES = if ("NES" %in% colnames(gsea_sig)) pw$NES else NA_real_
    )
  })

  if (nrow(tri_one) == 0) {
    skip_log[[folder_name]] <- "no overlap between top_genes and leadingEdge"
    next
  }

  triplets_list[[folder_name]] <- tri_one
}

triplets <- bind_rows(triplets_list)

if (nrow(triplets) == 0) {
  stop("没有成功构造 KO_gene -> affected_gene -> pathway 三层数据。")
}

if (length(skip_log) > 0) {
  skip_df <- data.frame(
    folder = names(skip_log),
    reason = unlist(skip_log),
    row.names = NULL
  )
  fwrite(skip_df, file.path(out_dir, "KO_skipped_reason.csv"))
}

shared_gene_summary <- triplets %>%
  distinct(KO_gene, affected_gene) %>%
  count(affected_gene, name = "n_KOs") %>%
  arrange(desc(n_KOs), affected_gene)

if (use_shared_gene_filter) {
  keep_genes <- shared_gene_summary %>%
    filter(n_KOs >= min_shared_kos) %>%
    pull(affected_gene)

  triplets <- triplets %>%
    filter(affected_gene %in% keep_genes)
}

if (nrow(triplets) == 0) {
  stop("共享基因过滤后没有剩余数据。")
}

top_pathways_keep <- triplets %>%
  group_by(pathway) %>%
  summarise(
    total_w = sum(gene_absZ),
    Count = n_distinct(affected_gene),
    minFDR = min(pathway_FDR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Count), desc(total_w), minFDR) %>%
  slice_head(n = show_top_n_pathways) %>%
  pull(pathway)

triplets <- triplets %>%
  filter(pathway %in% top_pathways_keep)

edges_ko_gene_all <- triplets %>%
  distinct(KO_gene, affected_gene, .keep_all = TRUE) %>%
  transmute(
    source = KO_gene,
    target = affected_gene,
    edge_type = "KO_to_gene",
    weight = gene_absZ
  )

edges_gene_path_all <- triplets %>%
  distinct(affected_gene, pathway, .keep_all = TRUE) %>%
  transmute(
    source = affected_gene,
    target = pathway,
    edge_type = "gene_to_pathway",
    weight = -log10(pathway_FDR + 1e-300)
  )

edges_all <- bind_rows(edges_ko_gene_all, edges_gene_path_all) %>%
  distinct()

ko_nodes <- triplets %>%
  distinct(KO_gene) %>%
  transmute(id = KO_gene, label = KO_gene, node_type = "KO_gene")

gene_nodes <- triplets %>%
  distinct(affected_gene) %>%
  transmute(id = affected_gene, label = affected_gene, node_type = "affected_gene")

pathway_nodes <- triplets %>%
  distinct(pathway) %>%
  transmute(id = pathway, label = pathway, node_type = "pathway")

nodes_all <- bind_rows(ko_nodes, gene_nodes, pathway_nodes) %>%
  distinct()

fwrite(edges_all, file.path(out_dir, "Tripartite_edges_for_Cytoscape.csv"))
fwrite(nodes_all, file.path(out_dir, "Tripartite_nodes_for_Cytoscape.csv"))

sankey_df <- triplets %>%
  distinct(KO_gene, affected_gene, pathway, .keep_all = TRUE) %>%
  mutate(
    axis1 = KO_gene,
    axis2 = affected_gene,
    axis3 = pathway,
    weight = gene_absZ
  )

ko_order <- sankey_df %>%
  group_by(KO_gene) %>%
  summarise(total_w = sum(weight), .groups = "drop") %>%
  arrange(desc(total_w)) %>%
  pull(KO_gene)

gene_order <- sankey_df %>%
  group_by(affected_gene) %>%
  summarise(total_w = sum(weight), .groups = "drop") %>%
  arrange(desc(total_w)) %>%
  pull(affected_gene)

pathway_order <- sankey_df %>%
  group_by(pathway) %>%
  summarise(
    total_w = sum(weight),
    Count = n_distinct(affected_gene),
    minFDR = min(pathway_FDR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(Count), desc(total_w), minFDR) %>%
  pull(pathway)

sankey_df <- sankey_df %>%
  mutate(
    axis1 = factor(axis1, levels = ko_order),
    axis2 = factor(axis2, levels = gene_order),
    axis3 = factor(axis3, levels = pathway_order)
  )

ko_palette   <- make_grad_desc("#0b3c49", length(ko_order),   light_col = "#eaf5f7")
gene_palette <- make_grad_desc("#2f6fb0", length(gene_order), light_col = "#edf4fb")
path_palette <- make_grad_desc("#6ec5e9", length(pathway_order), light_col = "#f2fbfe")

names(ko_palette)   <- ko_order
names(gene_palette) <- gene_order
names(path_palette) <- pathway_order

all_fill_values <- c(ko_palette, gene_palette, path_palette)

top_gene_labels <- sankey_df %>%
  distinct(axis2, .keep_all = TRUE) %>%
  arrange(desc(weight)) %>%
  slice_head(n = max_gene_labels_total) %>%
  pull(axis2) %>%
  as.character()

p_tmp <- ggplot(
  sankey_df,
  aes(
    y = weight,
    axis1 = axis1,
    axis2 = axis2,
    axis3 = axis3
  )
) +
  geom_alluvium(alpha = 0.2, width = 0.14) +
  geom_stratum(width = 0.14)

gb <- ggplot_build(p_tmp)
stratum_data <- gb$data[[2]]

stratum_left <- stratum_data %>%
  dplyr::filter(x == 1) %>%
  mutate(label = as.character(stratum))

stratum_mid <- stratum_data %>%
  dplyr::filter(x == 2) %>%
  mutate(label = as.character(stratum))

stratum_right <- stratum_data %>%
  dplyr::filter(x == 3) %>%
  mutate(label = as.character(stratum))

stratum_mid_show <- stratum_mid %>%
  filter(label %in% top_gene_labels)

y_min_all <- min(stratum_data$ymin, na.rm = TRUE)
y_max_all <- max(stratum_data$ymax, na.rm = TRUE)

n_ko <- length(ko_order)
total_height <- y_max_all - y_min_all
gap_frac <- 0.08
gap_abs <- total_height * gap_frac / n_ko
box_h <- (total_height - gap_abs * (n_ko - 1)) / n_ko

ko_equal_rect <- tibble(
  label = ko_order,
  x = 1,
  ymax = y_max_all - (seq_len(n_ko) - 1) * (box_h + gap_abs),
  ymin = y_max_all - (seq_len(n_ko)) * box_h - (seq_len(n_ko) - 1) * gap_abs
)

p_sankey <- ggplot(
  sankey_df,
  aes(
    y = weight,
    axis1 = axis1,
    axis2 = axis2,
    axis3 = axis3
  )
) +
  geom_alluvium(
    aes(fill = axis1),
    alpha = 0.28,
    width = 0.10,
    knot.pos = 0.35,
    color = NA
  ) +
  scale_x_discrete(
    limits = c("KO gene", "Affected gene", "Pathway"),
    expand = c(0.08, 0.08)
  ) +
  coord_cartesian(
    ylim = c(y_min_all, y_max_all),
    clip = "off"
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.margin = margin(0, 4, 0, 0)
  )

p_sankey <- p_sankey +
  geom_rect(
    data = ko_equal_rect,
    aes(
      xmin = x - 0.14, xmax = x + 0.14,
      ymin = ymin, ymax = ymax,
      fill = label
    ),
    inherit.aes = FALSE,
    color = "white",
    linewidth = 0.35
  )

p_sankey <- p_sankey +
  geom_rect(
    data = stratum_mid_show,
    aes(
      xmin = x - 0.14, xmax = x + 0.14,
      ymin = ymin, ymax = ymax,
      fill = label
    ),
    inherit.aes = FALSE,
    color = "white",
    linewidth = 0.35
  )

p_sankey <- p_sankey +
  geom_rect(
    data = stratum_right,
    aes(
      xmin = x - 0.14, xmax = x + 0.14,
      ymin = ymin, ymax = ymax,
      fill = label
    ),
    inherit.aes = FALSE,
    color = "white",
    linewidth = 0.35
  ) +
  scale_fill_manual(
    values = all_fill_values,
    guide = "none"
  )

p_sankey <- p_sankey +
  geom_text(
    data = ko_equal_rect,
    aes(
      x = x - 0.16,
      y = (ymin + ymax) / 2,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1,
    size = ko_label_size
  ) +
  geom_text(
    data = stratum_mid_show,
    aes(
      x = x - 0.16,
      y = (ymin + ymax) / 2,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1,
    size = gene_label_size
  ) +
  geom_text(
    data = stratum_right,
    aes(
      x = x - 0.16,
      y = (ymin + ymax) / 2,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1,
    size = pathway_label_size
  )

bubble_df <- sankey_df %>%
  group_by(axis3) %>%
  summarise(
    Count = n_distinct(axis2),
    pvalue = min(pathway_FDR, na.rm = TRUE),
    total_w = sum(weight),
    .groups = "drop"
  ) %>%
  rename(pathway = axis3) %>%
  mutate(pathway = as.character(pathway))

right_pos <- stratum_right %>%
  transmute(
    pathway = as.character(label),
    y_center = (ymin + ymax) / 2
  )

bubble_df <- bubble_df %>%
  inner_join(right_pos, by = "pathway") %>%
  mutate(
    xval = -log10(pvalue + 1e-300)
  ) %>%
  arrange(desc(y_center))

xmax_bubble <- max(bubble_df$xval, na.rm = TRUE)

p_bubble <- ggplot(bubble_df, aes(x = xval, y = y_center)) +
  geom_point(
    aes(size = Count, fill = xval),
    shape = 21,
    color = "#1a1a1a",
    stroke = 0.25,
    alpha = 0.95
  ) +
  scale_size_continuous(
    range = c(3.5, 10),
    name = "Gene count"
  ) +
  scale_fill_gradient(
    low = "#dff4fb",
    high = "#6ec5e9",
    name = "-log10(FDR)"
  ) +
  scale_x_continuous(
    limits = c(0, xmax_bubble * 1.10),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  coord_cartesian(
    ylim = c(y_min_all, y_max_all),
    clip = "off"
  ) +
  labs(
    x = "-log10(FDR)",
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.box = "vertical",
    plot.margin = margin(0, 10, 0, 0)
  )

legend_grob <- cowplot::get_legend(p_bubble)
p_bubble_noleg <- p_bubble + theme(legend.position = "none")

main_body <- cowplot::plot_grid(
  p_sankey,
  p_bubble_noleg,
  nrow = 1,
  rel_widths = c(3.45, 1.08),
  align = "v",
  axis = "lr"
)

right_col <- cowplot::plot_grid(
  legend_grob,
  NULL,
  ncol = 1,
  rel_heights = c(1.0, 2.8)
)

p_final <- cowplot::plot_grid(
  main_body,
  right_col,
  nrow = 1,
  rel_widths = c(5.25, 0.95)
)

ggsave(
  filename = file.path(out_dir, "KO_gene_affected_gene_pathway_sankey_bubble_v9.pdf"),
  plot = p_final,
  width = pdf_width,
  height = pdf_height
)

ggsave(
  filename = file.path(out_dir, "KO_gene_affected_gene_pathway_sankey_bubble_v9.png"),
  plot = p_final,
  width = pdf_width,
  height = pdf_height,
  dpi = 300
)

fwrite(triplets, file.path(out_dir, "sankey_triplets_used.csv"))

cat(
  "输出完成：\n",
  file.path(out_dir, "KO_gene_affected_gene_pathway_sankey_bubble_v9.pdf"), "\n",
  file.path(out_dir, "KO_gene_affected_gene_pathway_sankey_bubble_v9.png"), "\n",
  file.path(out_dir, "Tripartite_edges_for_Cytoscape.csv"), "\n",
  file.path(out_dir, "Tripartite_nodes_for_Cytoscape.csv"), "\n",
  file.path(out_dir, "sankey_triplets_used.csv"), "\n"
)
