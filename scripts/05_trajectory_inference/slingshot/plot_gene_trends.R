rm(list = ls()); gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(data.table)
})

rds_path <- project_path("results_root", "single_cell/slingshot/rds/preHTC_HTC_Seurat_with_slingshot.rds")
stopifnot(file.exists(rds_path))

outdir <- project_path("results_root", "single_cell/slingshot/pseudotime_geneplots_L1")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

single_dir <- file.path(outdir, "single_gene_plots")
dir.create(single_dir, showWarnings = FALSE, recursive = TRUE)

pseudotime_col <- "pseudotime_L1"

genes_use <- c(
  "MDH1","DCN","PKM","LDHA","SLC35A3","PLOD2","CD44","ENO1","PPIA","AK3",
  "TGFBI","AGL","NASP","LHX9","COPB2","PAM","VCAN","ELF3","SDC2","ME2",
  "PGK1","HDLBP","COG2"
)

bin_width <- 1
min_cells_per_bin <- 5

restrict_celltype <- FALSE
celltype_col <- "cell_type"
celltypes_keep <- c("preHTC", "HTC")

message("Reading Seurat object: ", rds_path)
seurat_obj <- readRDS(rds_path)
stopifnot(inherits(seurat_obj, "Seurat"))

if (!pseudotime_col %in% colnames(seurat_obj@meta.data)) {
  stop("对象中不存在列: ", pseudotime_col,
       "\n可用 pseudotime 列为: ",
       paste(grep("^pseudotime", colnames(seurat_obj@meta.data), value = TRUE), collapse = ", "))
}

if (restrict_celltype) {
  if (!celltype_col %in% colnames(seurat_obj@meta.data)) {
    stop("meta.data 中不存在列: ", celltype_col)
  }
  keep_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[celltype_col]] %in% celltypes_keep]
  if (length(keep_cells) == 0) {
    stop("未找到目标 cell_type: ", paste(celltypes_keep, collapse = ", "))
  }
  seurat_obj <- subset(seurat_obj, cells = keep_cells)
  message("Subset done. Cells kept: ", ncol(seurat_obj))
}

all_features <- rownames(seurat_obj)
genes_present <- intersect(genes_use, all_features)
genes_missing <- setdiff(genes_use, all_features)

writeLines(genes_present, file.path(outdir, "genes_present.txt"))
writeLines(genes_missing, file.path(outdir, "genes_missing.txt"))

if (length(genes_present) == 0) {
  stop("目标基因均未在对象中找到。")
}

message("Genes found: ", length(genes_present))
if (length(genes_missing) > 0) {
  message("Missing genes: ", paste(genes_missing, collapse = ", "))
}

vars_use <- c(genes_present, pseudotime_col)
raw_data <- FetchData(seurat_obj, vars = vars_use)
colnames(raw_data)[ncol(raw_data)] <- "Pseudotime"

raw_data <- raw_data %>%
  mutate(cell = rownames(.)) %>%
  filter(!is.na(Pseudotime))

long_data <- raw_data %>%
  pivot_longer(cols = all_of(genes_present), names_to = "Gene", values_to = "Expression")

fwrite(as.data.table(long_data), file.path(outdir, "raw_long_expression_with_pseudotime.csv"))

plot_data_binned <- long_data %>%
  mutate(Bin_Start = floor(Pseudotime / bin_width) * bin_width) %>%
  mutate(Bin_Label = paste0(Bin_Start, "-", Bin_Start + bin_width)) %>%
  group_by(Gene, Bin_Start, Bin_Label) %>%
  summarise(
    Mean_Expression = mean(Expression, na.rm = TRUE),
    Cell_Count = n(),
    .groups = "drop"
  ) %>%
  filter(Cell_Count >= min_cells_per_bin) %>%
  arrange(Gene, Bin_Start)

fwrite(as.data.table(plot_data_binned), file.path(outdir, "binned_mean_expression.csv"))

plot_data_final <- plot_data_binned %>%
  group_by(Gene) %>%
  mutate(
    Scaled_Expression = rescale(Mean_Expression, to = c(0, 1))
  ) %>%
  ungroup()

bin_levels <- plot_data_final %>%
  distinct(Bin_Start, Bin_Label) %>%
  arrange(Bin_Start) %>%
  pull(Bin_Label)

plot_data_final$Bin_Label <- factor(plot_data_final$Bin_Label, levels = bin_levels)

all_labels <- levels(plot_data_final$Bin_Label)
if (length(all_labels) <= 12) {
  breaks_to_show <- all_labels
} else {
  breaks_to_show <- all_labels[seq(1, length(all_labels), by = 3)]
}

fwrite(as.data.table(plot_data_final), file.path(outdir, "plot_data_final_scaled.csv"))

p_all <- ggplot(
  plot_data_final,
  aes(x = Bin_Label, y = Scaled_Expression, color = Gene, group = Gene)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.8) +
  theme_classic(base_size = 14) +
  scale_x_discrete(breaks = breaks_to_show) +
  labs(
    title = "Dynamic changes of candidate genes along pseudotime",
    subtitle = paste0("Scaled mean expression by pseudotime bins (", pseudotime_col, ")"),
    x = "Pseudotime bins",
    y = "Relative expression (scaled 0-1)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "right",
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
  )

ggsave(
  filename = file.path(outdir, "All_genes_pseudotime_L1_overview.pdf"),
  plot = p_all, width = 12, height = 7
)
ggsave(
  filename = file.path(outdir, "All_genes_pseudotime_L1_overview.png"),
  plot = p_all, width = 12, height = 7, dpi = 300
)

p_facet <- ggplot(
  plot_data_final,
  aes(x = Bin_Label, y = Scaled_Expression, group = 1)
) +
  geom_col(fill = "#9ecae1", alpha = 0.75, width = 0.75) +
  geom_line(color = "#08519c", linewidth = 0.9) +
  geom_point(color = "#08306b", size = 1.8) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  scale_x_discrete(breaks = breaks_to_show) +
  theme_classic(base_size = 13) +
  labs(
    title = "Expression dynamics of candidate genes along pseudotime",
    subtitle = paste0("Binned mean expression using ", pseudotime_col),
    x = "Pseudotime bins",
    y = "Relative expression (scaled 0-1)"
  ) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 9),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
  )

ggsave(
  filename = file.path(outdir, "Facet_genes_pseudotime_L1.pdf"),
  plot = p_facet, width = 16, height = 12
)
ggsave(
  filename = file.path(outdir, "Facet_genes_pseudotime_L1.png"),
  plot = p_facet, width = 16, height = 12, dpi = 300
)

for (g in unique(plot_data_final$Gene)) {

  df_g <- plot_data_final %>% filter(Gene == g)

  p_g <- ggplot(df_g, aes(x = Bin_Label, y = Scaled_Expression, group = 1)) +
    geom_col(fill = "#4e79a7", alpha = 0.70, width = 0.7) +
    geom_line(color = "#1f3b5d", linewidth = 1.2) +
    geom_point(size = 2.2, shape = 21, fill = "white", color = "#1f3b5d", stroke = 1.2) +
    scale_x_discrete(breaks = breaks_to_show) +
    theme_classic(base_size = 15) +
    labs(
      title = paste0(g, " along pseudotime"),
      subtitle = paste0("Scaled mean expression (", pseudotime_col, ")"),
      x = "Pseudotime bins",
      y = "Relative expression (scaled 0-1)"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
    )

  ggsave(
    filename = file.path(single_dir, paste0(g, "_pseudotime_L1.pdf")),
    plot = p_g, width = 8, height = 5.5
  )
  ggsave(
    filename = file.path(single_dir, paste0(g, "_pseudotime_L1.png")),
    plot = p_g, width = 8, height = 5.5, dpi = 300
  )
}

safe_spearman_p <- function(x, y) {
  out <- tryCatch(
    suppressWarnings(cor.test(x, y, method = "spearman")$p.value),
    error = function(e) NA_real_
  )
  out
}

trend_tbl <- long_data %>%
  group_by(Gene) %>%
  summarise(
    Spearman_rho = suppressWarnings(cor(Expression, Pseudotime, method = "spearman", use = "complete.obs")),
    Spearman_p = safe_spearman_p(Expression, Pseudotime),
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(Spearman_p, method = "BH")) %>%
  arrange(desc(Spearman_rho))

fwrite(as.data.table(trend_tbl), file.path(outdir, "gene_pseudotime_spearman_summary.csv"))

print(p_all)
print(p_facet)

message("Done.")
message("Input RDS: ", rds_path)
message("Output dir: ", outdir)
message("Pseudotime column: ", pseudotime_col)
message("Present genes: ", paste(genes_present, collapse = ", "))
if (length(genes_missing) > 0) {
  message("Missing genes: ", paste(genes_missing, collapse = ", "))
}
