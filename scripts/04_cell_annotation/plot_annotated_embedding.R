library(ggpubr)
library(dplyr)

plot_data <- cbind(
  seurat_object@reductions$umap@cell.embeddings,
  Subcluster = seurat_object$seurat_clusters
) %>%
  as.data.frame()

colnames(plot_data)[1:2] <- c("UMAP_1", "UMAP_2")

plot_data$Subcluster <- paste0("P", plot_data$Subcluster)

cols1 <- c(
  "#bae4c6",
  "#72a1df",
  "#ead5f4",
  "#baa28b",
  "#f9f897",
  "#64dec8",
  "#a78dd3",
  "#d48073"

)

p1 <- ggscatter(
  plot_data,
  x        = "UMAP_1",
  y        = "UMAP_2",
  color    = "Subcluster",
  alpha    = 1,
  size     = 1,
  ellipse  = TRUE,
  palette  = rev(cols1),
  ggtheme  = theme_classic2() +
    theme(
      text = element_text(
        family = "serif",
        face   = "bold",
        size   = 13
      )
    )
)

p1
