cran_packages <- c(
  "BisqueRNA", "cowplot", "data.table", "dplyr", "ggalluvial",
  "ggplot2", "ggpubr", "harmony", "Matrix", "matrixStats",
  "msigdbr", "Nebulosa", "patchwork", "pheatmap", "pROC",
  "purrr", "RColorBrewer", "reshape2", "rlang", "rTensor",
  "scales", "Seurat", "SeuratObject", "stringr", "tibble",
  "tidyr", "tidyverse", "uwot", "VennDiagram", "viridis",
  "wesanderson"
)

bioconductor_packages <- c(
  "AnnotationDbi", "AUCell", "Biobase", "BiocGenerics",
  "BiocParallel", "celldex", "clusterProfiler", "DESeq2", "edgeR",
  "enrichplot", "fgsea", "limma", "org.Hs.eg.db", "org.Mm.eg.db",
  "S4Vectors", "SingleCellExperiment", "SingleR", "slingshot",
  "SummarizedExperiment", "tradeSeq", "UCell"
)

specialized_packages <- c(
  "CytoTRACE2", "DoubletFinder", "monocle3", "MuSiC", "SCopeLoomR",
  "SCP", "scrapper", "scTenifoldKnk", "SeuratWrappers"
)

install_oagly_dependencies <- function(install_specialized = FALSE) {
  missing_cran <- setdiff(cran_packages, rownames(installed.packages()))
  if (length(missing_cran)) {
    install.packages(missing_cran, repos = "https://cloud.r-project.org")
  }

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  missing_bioc <- setdiff(bioconductor_packages, rownames(installed.packages()))
  if (length(missing_bioc)) {
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  }

  if (install_specialized) {
    message(
      "Install specialized packages according to their upstream documentation: ",
      paste(specialized_packages, collapse = ", ")
    )
  }
}
