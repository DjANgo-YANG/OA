library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)

meta <- seurat_object@meta.data

cond_col <- "condition"
x_col    <- "region"
fill_col <- "cell_type2"

stopifnot(cond_col %in% colnames(meta),
          x_col    %in% colnames(meta),
          fill_col %in% colnames(meta))

meta[[x_col]] <- factor(meta[[x_col]], levels = c("AS","SZ","MZ","DZ"))

fill_levels <- meta %>% pull(.data[[fill_col]]) %>% as.character() %>% unique() %>% sort()
meta[[fill_col]] <- factor(meta[[fill_col]], levels = fill_levels)

plot_stacked_prop <- function(meta_sub, title_text) {
  df <- meta_sub %>%
    count(.data[[x_col]], .data[[fill_col]], name = "n") %>%
    group_by(.data[[x_col]]) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()

  ggplot(df, aes(x = .data[[x_col]], y = prop, fill = .data[[fill_col]])) +
    geom_col(width = 0.85) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(title = title_text, x = x_col, y = "Proportion", fill = fill_col) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
}

plot_stacked_count <- function(meta_sub, title_text) {
  df <- meta_sub %>% count(.data[[x_col]], .data[[fill_col]], name = "n")
  ggplot(df, aes(x = .data[[x_col]], y = n, fill = .data[[fill_col]])) +
    geom_col(width = 0.85) +
    labs(title = title_text, x = x_col, y = "Cell count", fill = fill_col) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
}

meta_HC <- meta %>% filter(.data[[cond_col]] == "HC")
meta_OA <- meta %>% filter(.data[[cond_col]] == "OA")

p_HC <- plot_stacked_prop(meta_HC, "HC")
p_OA <- plot_stacked_prop(meta_OA, "OA")

p_all <- (p_HC + p_OA) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

p_all

ggsave("stacked_HC_vs_OA.pdf", p_all, width = 10, height = 4)
