library(dplyr)
library(stringr)

props_long <- props_long %>%
  mutate(

    group = case_when(
      str_starts(sample, "Control") ~ "HC",
      str_starts(sample, "NWB")     ~ "OA",
      str_starts(sample, "WB")      ~ "OA",
      TRUE                          ~ NA_character_
    ),

    region = case_when(
      str_detect(sample, "\\.AS") ~ "AS",
      str_detect(sample, "\\.SZ") ~ "SZ",
      str_detect(sample, "\\.MZ") ~ "MZ",
      str_detect(sample, "\\.DZ") ~ "DZ",
      TRUE                        ~ NA_character_
    )
  )

props_long <- props_long %>%
  filter(!is.na(group), !is.na(region))

celltype_levels <- c("EC", "HomC", "HTC",
                     "MMP13-", "MMP13+",
                     "preHTC", "proC",
                     "RegC", "RepC")

props_long <- props_long %>%
  mutate(
    cell_type = factor(cell_type, levels = celltype_levels),
    group     = factor(group,     levels = c("HC", "OA")),
    region    = factor(region,    levels = c("AS", "SZ", "MZ", "DZ"))
  )
table(props_long$group, props_long$region)
library(ggplot2)
library(scales)

props_group <- props_long %>%
  group_by(group, cell_type) %>%
  summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop")

p_group_rel <- ggplot(props_group,
                      aes(x = group, y = mean_prop, fill = cell_type)) +
  geom_col(position = "fill", color = NA) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Group",
    y = "Relative proportion of cell states",
    fill = "Cell state"
  ) +
  theme_bw() +
  theme(
    panel.grid  = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

p_group_rel
props_region <- props_long %>%
  group_by(group, region, cell_type) %>%
  summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop")

p_region_rel <- ggplot(props_region,
                       aes(x = region, y = mean_prop, fill = cell_type)) +
  geom_col(position = "fill", color = NA) +
  facet_grid(. ~ group) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Region",
    y = "Relative proportion of cell states",
    fill = "Cell state"
  ) +
  theme_bw() +
  theme(
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text       = element_text(size = 12, face = "bold"),
    axis.text.x      = element_text(angle = 0, hjust = 0.5)
  )

p_region_rel
p_region_abs <- ggplot(props_region,
                       aes(x = region, y = mean_prop, fill = cell_type)) +
  geom_col(position = "stack", color = NA) +
  facet_grid(. ~ group) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Region",
    y = "Average fraction (absolute)",
    fill = "Cell state"
  ) +
  theme_bw() +
  theme(
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "grey90", colour = "black"),
    strip.text       = element_text(size = 12, face = "bold"),
    axis.text.x      = element_text(angle = 0, hjust = 0.5)
  )

p_region_abs
