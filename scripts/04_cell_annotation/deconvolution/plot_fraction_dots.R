source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)
library(scales)

prop_file <- project_path("data_root", "bulk_deconvolution/Bisque_cell_fraction_matrix.csv")

prop_mat <- read.csv(prop_file, row.names = 1, check.names = FALSE)
prop_mat <- as.matrix(prop_mat)

cat("Top rownames:\n"); print(head(rownames(prop_mat), 10))
cat("Top colnames:\n"); print(head(colnames(prop_mat), 10))

celltype_like <- function(x) any(x %in% c("EC","HomC","HTC","MMP13-","MMP13+","preHTC","proC","RegC","RepC"))
sample_like   <- function(x) any(grepl("(Control|HC|OA|WB|NWB|AS|SZ|MZ|DZ)", x, ignore.case = TRUE))

need_transpose <- celltype_like(rownames(prop_mat)) && sample_like(colnames(prop_mat))

if (need_transpose) {
  cat(">>> Detected: rows look like cell types and cols look like samples. Transposing matrix...\n")
  prop_mat <- t(prop_mat)
}

cat("After fix, example sample names:\n")
print(head(rownames(prop_mat), 20))
cat("After fix, cell types:\n")
print(colnames(prop_mat))

props_long <- as.data.frame(prop_mat) %>%
  rownames_to_column("sample") %>%
  pivot_longer(cols = -sample, names_to = "cell_type", values_to = "prop") %>%
  mutate(
    region = str_extract(sample, "(AS|SZ|MZ|DZ)"),
    group = case_when(
      str_detect(sample, regex("Control|HC", ignore_case = TRUE)) ~ "HC",
      str_detect(sample, regex("WB|NWB|OA", ignore_case = TRUE)) ~ "OA",
      TRUE ~ NA_character_
    )
  )

bad_samples <- props_long %>%
  distinct(sample, group, region) %>%
  filter(is.na(group) | is.na(region))

cat("Bad samples (cannot parse group/region): ", nrow(bad_samples), "\n")
if (nrow(bad_samples) > 0) print(head(bad_samples, 50))

props_long <- props_long %>%
  filter(!is.na(group), !is.na(region)) %>%
  mutate(
    group  = factor(group,  levels = c("HC", "OA")),
    region = factor(region, levels = c("AS", "SZ", "MZ", "DZ"))
  )

stopifnot(nrow(props_long) > 0)

mean_df <- props_long %>%
  group_by(group, region, cell_type) %>%
  summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop")

stat_df <- props_long %>%
  group_by(region, cell_type) %>%
  summarise(
    mean_HC = mean(prop[group == "HC"], na.rm = TRUE),
    mean_OA = mean(prop[group == "OA"], na.rm = TRUE),
    p = {
      x <- prop[group == "HC"]; y <- prop[group == "OA"]
      if (length(na.omit(x)) < 2 || length(na.omit(y)) < 2) NA_real_
      else suppressWarnings(wilcox.test(x, y)$p.value)
    },
    .groups = "drop"
  ) %>%
  mutate(
    diff = mean_OA - mean_HC,
    signed_logp = ifelse(is.na(p), 0, sign(diff) * pmin(3, -log10(p)))
  )

plot_df <- mean_df %>%
  left_join(stat_df %>% select(region, cell_type, signed_logp),
            by = c("region", "cell_type"))

p <- ggplot(plot_df, aes(x = cell_type, y = region)) +
  geom_point(aes(size = mean_prop, color = signed_logp)) +
  facet_grid(group ~ ., switch = "y") +
  scale_size(range = c(0.5, 8), labels = percent_format(accuracy = 1)) +
  scale_color_gradient2(
    low = "#9769df", mid = "#69c9df", high = "#df7f69",
    limits = c(-3, 3), oob = scales::squish,
    name = expression(paste("signed  ", -log[10], "(p)"))
  ) +
  labs(x = NULL, y = NULL, size = "Mean proportion") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey92", colour = "black"),
    strip.text = element_text(face = "bold")
  )

print(p)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)
library(scales)

props_long <- as.data.frame(prop_mat) %>%
  rownames_to_column("sample") %>%
  pivot_longer(cols = -sample, names_to = "cell_type", values_to = "prop") %>%
  mutate(
    group = case_when(
      str_detect(sample, regex("Control|HC", ignore_case = TRUE)) ~ "HC",
      str_detect(sample, regex("WB|NWB|OA", ignore_case = TRUE)) ~ "OA",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = c("HC", "OA"))
  ) %>%
  filter(!is.na(group))

stopifnot(nrow(props_long) > 0)

stat_df <- props_long %>%
  group_by(cell_type) %>%
  summarise(
    mean_HC = mean(prop[group == "HC"], na.rm = TRUE),
    mean_OA = mean(prop[group == "OA"], na.rm = TRUE),
    p = {
      x <- prop[group == "HC"]; y <- prop[group == "OA"]
      if (length(na.omit(x)) < 2 || length(na.omit(y)) < 2) NA_real_
      else suppressWarnings(wilcox.test(x, y)$p.value)
    },
    .groups = "drop"
  ) %>%
  mutate(
    diff = mean_OA - mean_HC,
    signed_logp = ifelse(is.na(p), 0, sign(diff) * pmin(3, -log10(p)))
  )

plot_df <- props_long %>%
  group_by(group, cell_type) %>%
  summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop") %>%
  left_join(stat_df %>% select(cell_type, signed_logp), by = "cell_type")

p <- ggplot(plot_df, aes(x = cell_type, y = group)) +
  geom_point(aes(size = mean_prop, color = signed_logp)) +
  scale_size(range = c(0.5, 8), labels = percent_format(accuracy = 1)) +
  scale_color_gradient2(
    low = "#341e45", mid = "#2d8f84", high = "#ecdc33",
    limits = c(-3, 3), oob = squish,
    name = expression(paste("signed  ", -log[10], "(p)"))
  ) +
  labs(x = NULL, y = NULL, size = "Mean proportion") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold")
  )

print(p)
