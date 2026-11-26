library(dplyr)
library(ggplot2)
library(scales)

obj <- seurat_object

group_col   <- "group"
subtype_col <- "cell_type"

meta <- obj@meta.data

meta_use <- meta %>%
  filter(!is.na(.data[[group_col]]),
         !is.na(.data[[subtype_col]]))

meta_use[[group_col]] <- factor(meta_use[[group_col]],
                                levels = c("HC", "OA"))

tab_cnt <- table(meta_use[[group_col]], meta_use[[subtype_col]])

df_cnt <- as.data.frame(tab_cnt)
colnames(df_cnt) <- c("group", "cell_subtype", "n")

p_abs <- ggplot(df_cnt, aes(x = group, y = n, fill = cell_subtype)) +
  geom_col(position = "stack") +
  ylab("Number of cells") +
  xlab("Group") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

df_prop <- df_cnt %>%
  group_by(group) %>%
  mutate(prop = n / sum(n))

p_rel <- ggplot(df_prop, aes(x = group, y = prop, fill = cell_subtype)) +
  geom_col(position = "fill") +
  ylab("Proportion of cells") +
  xlab("Group") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

p_abs
p_rel
