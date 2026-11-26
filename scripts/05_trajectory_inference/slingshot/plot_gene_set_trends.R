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

gmt_inflam <- project_path("references", "gene_sets/HALLMARK_INFLAMMATORY_RESPONSE.v2026.1.Hs.gmt")
gmt_ecmdeg <- project_path("references", "gene_sets/REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX.v2025.1.Hs.gmt")

stopifnot(file.exists(gmt_inflam))
stopifnot(file.exists(gmt_ecmdeg))

outdir <- project_path("results_root", "single_cell/slingshot/pseudotime_genesetplots_L1")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

pseudotime_col <- "pseudotime_L1"
bin_width <- 2
min_cells_per_bin <- 5

genes_23 <- c(
  "MDH1","DCN","PKM","LDHA","SLC35A3","PLOD2","CD44","ENO1","PPIA","AK3",
  "TGFBI","AGL","NASP","LHX9","COPB2","PAM","VCAN","ELF3","SDC2","ME2",
  "PGK1","HDLBP","COG2"
)

read_gmt_simple <- function(gmt_file) {
  x <- readLines(gmt_file, warn = FALSE, encoding = "UTF-8")
  x <- x[nzchar(x)]
  spl <- strsplit(x, "\t", fixed = TRUE)
  gs_list <- lapply(spl, function(z) unique(z[-c(1,2)]))
  gs_name <- vapply(spl, `[`, character(1), 1)
  names(gs_list) <- gs_name
  gs_list
}

safe_spearman_p <- function(x, y) {
  tryCatch(
    suppressWarnings(cor.test(x, y, method = "spearman")$p.value),
    error = function(e) NA_real_
  )
}

message("Reading Seurat object: ", rds_path)
seurat_obj <- readRDS(rds_path)
stopifnot(inherits(seurat_obj, "Seurat"))

if (!pseudotime_col %in% colnames(seurat_obj@meta.data)) {
  stop("对象中不存在列: ", pseudotime_col)
}

gs_inflam_raw <- read_gmt_simple(gmt_inflam)
gs_ecmdeg_raw <- read_gmt_simple(gmt_ecmdeg)

gs_inflam_name <- names(gs_inflam_raw)[1]
gs_ecmdeg_name <- names(gs_ecmdeg_raw)[1]

gs_inflam <- unique(gs_inflam_raw[[1]])
gs_ecmdeg <- unique(gs_ecmdeg_raw[[1]])

all_features <- rownames(seurat_obj)

gs_inflam_use <- intersect(gs_inflam, all_features)
gs_ecmdeg_use <- intersect(gs_ecmdeg, all_features)
gs_23_use     <- intersect(genes_23, all_features)

if (length(gs_inflam_use) < 5) {
  stop("HALLMARK_INFLAMMATORY_RESPONSE 在对象中匹配到的基因太少: ", length(gs_inflam_use))
}
if (length(gs_ecmdeg_use) < 5) {
  stop("REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX 在对象中匹配到的基因太少: ", length(gs_ecmdeg_use))
}
if (length(gs_23_use) < 5) {
  stop("23基因集在对象中匹配到的基因太少: ", length(gs_23_use))
}

writeLines(gs_inflam_use, file.path(outdir, "geneset_inflammatory_present.txt"))
writeLines(gs_ecmdeg_use, file.path(outdir, "geneset_ecmdeg_present.txt"))
writeLines(gs_23_use,     file.path(outdir, "geneset_23_present.txt"))

writeLines(setdiff(gs_inflam, all_features), file.path(outdir, "geneset_inflammatory_missing.txt"))
writeLines(setdiff(gs_ecmdeg, all_features), file.path(outdir, "geneset_ecmdeg_missing.txt"))
writeLines(setdiff(genes_23, all_features),  file.path(outdir, "geneset_23_missing.txt"))

set.seed(1234)
seurat_obj <- AddModuleScore(
  object = seurat_obj,
  features = list(gs_inflam_use),
  name = "InflammatoryScore"
)

set.seed(1234)
seurat_obj <- AddModuleScore(
  object = seurat_obj,
  features = list(gs_ecmdeg_use),
  name = "ECMDegradationScore"
)

set.seed(1234)
seurat_obj <- AddModuleScore(
  object = seurat_obj,
  features = list(gs_23_use),
  name = "Gene23Score"
)

score_cols <- c("InflammatoryScore1", "ECMDegradationScore1", "Gene23Score1")
score_labels <- c(
  "Inflammatory_Response",
  "ECM_Degradation",
  "GeneSet_23"
)

raw_data <- FetchData(seurat_obj, vars = c(score_cols, pseudotime_col))
colnames(raw_data) <- c(score_labels, "Pseudotime")
raw_data <- raw_data %>%
  mutate(cell = rownames(.)) %>%
  filter(!is.na(Pseudotime))

long_data <- raw_data %>%
  pivot_longer(
    cols = all_of(score_labels),
    names_to = "GeneSet",
    values_to = "Score"
  )

fwrite(as.data.table(long_data), file.path(outdir, "raw_geneset_scores_with_pseudotime.csv"))

plot_data_binned <- long_data %>%
  mutate(Bin_Start = floor(Pseudotime / bin_width) * bin_width) %>%
  mutate(Bin_Label = paste0(Bin_Start, "-", Bin_Start + bin_width)) %>%
  group_by(GeneSet, Bin_Start, Bin_Label) %>%
  summarise(
    Mean_Score = mean(Score, na.rm = TRUE),
    Cell_Count = n(),
    .groups = "drop"
  ) %>%
  filter(Cell_Count >= min_cells_per_bin) %>%
  arrange(GeneSet, Bin_Start)

fwrite(as.data.table(plot_data_binned), file.path(outdir, "binned_mean_geneset_scores.csv"))

plot_data_final <- plot_data_binned %>%
  group_by(GeneSet) %>%
  mutate(
    Scaled_Score = rescale(Mean_Score, to = c(0, 1))
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
  aes(x = Bin_Label, y = Scaled_Score, color = GeneSet, group = GeneSet)
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.2) +
  theme_classic(base_size = 14) +
  scale_x_discrete(breaks = breaks_to_show) +
  labs(
    title = "Gene set dynamics along pseudotime",
    subtitle = paste0("Binned mean module scores using ", pseudotime_col),
    x = "Pseudotime bins",
    y = "Relative score (scaled 0-1)",
    color = "Gene set"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "right",
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
  )

ggsave(
  filename = file.path(outdir, "GeneSets_pseudotime_L1_overview.pdf"),
  plot = p_all, width = 10, height = 6
)
ggsave(
  filename = file.path(outdir, "GeneSets_pseudotime_L1_overview.png"),
  plot = p_all, width = 10, height = 6, dpi = 300
)

p_facet <- ggplot(
  plot_data_final,
  aes(x = Bin_Label, y = Scaled_Score, group = 1)
) +
  geom_col(fill = "#9ecae1", alpha = 0.75, width = 0.75) +
  geom_line(color = "#08519c", linewidth = 1.0) +
  geom_point(color = "#08306b", size = 2) +
  facet_wrap(~ GeneSet, scales = "free_y", ncol = 1) +
  scale_x_discrete(breaks = breaks_to_show) +
  theme_classic(base_size = 14) +
  labs(
    title = "Gene set score dynamics along pseudotime",
    subtitle = paste0("Binned mean module scores using ", pseudotime_col),
    x = "Pseudotime bins",
    y = "Relative score (scaled 0-1)"
  ) +
  theme(
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
  )

ggsave(
  filename = file.path(outdir, "GeneSets_pseudotime_L1_facet.pdf"),
  plot = p_facet, width = 10, height = 10
)
ggsave(
  filename = file.path(outdir, "GeneSets_pseudotime_L1_facet.png"),
  plot = p_facet, width = 10, height = 10, dpi = 300
)

single_dir <- file.path(outdir, "single_geneset_plots")
dir.create(single_dir, showWarnings = FALSE, recursive = TRUE)

for (gs in unique(plot_data_final$GeneSet)) {

  df_gs <- plot_data_final %>% filter(GeneSet == gs)

  p_gs <- ggplot(df_gs, aes(x = Bin_Label, y = Scaled_Score, group = 1)) +
    geom_col(fill = "#4e79a7", alpha = 0.70, width = 0.7) +
    geom_line(color = "#1f3b5d", linewidth = 1.2) +
    geom_point(size = 2.3, shape = 21, fill = "white", color = "#1f3b5d", stroke = 1.2) +
    scale_x_discrete(breaks = breaks_to_show) +
    theme_classic(base_size = 15) +
    labs(
      title = paste0(gs, " along pseudotime"),
      subtitle = paste0("Scaled mean module score (", pseudotime_col, ")"),
      x = "Pseudotime bins",
      y = "Relative score (scaled 0-1)"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")
    )

  ggsave(
    filename = file.path(single_dir, paste0(gs, "_pseudotime_L1.pdf")),
    plot = p_gs, width = 8, height = 5.5
  )
  ggsave(
    filename = file.path(single_dir, paste0(gs, "_pseudotime_L1.png")),
    plot = p_gs, width = 8, height = 5.5, dpi = 300
  )
}

trend_tbl <- long_data %>%
  group_by(GeneSet) %>%
  summarise(
    Spearman_rho = suppressWarnings(cor(Score, Pseudotime, method = "spearman", use = "complete.obs")),
    Spearman_p = safe_spearman_p(Score, Pseudotime),
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(Spearman_p, method = "BH")) %>%
  arrange(desc(Spearman_rho))

fwrite(as.data.table(trend_tbl), file.path(outdir, "geneset_pseudotime_spearman_summary.csv"))

geneset_info <- data.frame(
  GeneSet = c("Inflammatory_Response", "ECM_Degradation", "GeneSet_23"),
  Source = c(gmt_inflam, gmt_ecmdeg, "Custom_23_genes"),
  N_genes_in_object = c(length(gs_inflam_use), length(gs_ecmdeg_use), length(gs_23_use)),
  stringsAsFactors = FALSE
)
write.csv(geneset_info, file.path(outdir, "geneset_info.csv"), row.names = FALSE)

print(p_all)
print(p_facet)

message("Done.")
message("Input RDS: ", rds_path)
message("Output dir: ", outdir)
message("Pseudotime column: ", pseudotime_col)
message("Inflammatory genes matched: ", length(gs_inflam_use))
message("ECM degradation genes matched: ", length(gs_ecmdeg_use))
message("23-gene set matched: ", length(gs_23_use))
