rm(list = ls())
gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))

library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(DoubletFinder)
library(future)
library(dplyr)

work_dir <- project_path("results_root", "single_cell/doublets")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)

seurat_object <- readRDS(project_path("results_root", "single_cell/merged/1.rds"))

future::plan("sequential")

seurat_list <- SplitObject(seurat_object, split.by = "orig.ident")

rm(seurat_object)
gc()

run_doubletfinder_one <- function(obj,
                                  PCs_use   = 1:30,
                                  nfeatures = 3000,
                                  exp_rate  = 0.07,
                                  pN        = 0.25)
{
  sample_id <- unique(obj$orig.ident)
  message("===== Running DoubletFinder on sample: ", sample_id, " =====")
  message("细胞数：", ncol(obj))

  DefaultAssay(obj) <- "RNA"

  obj <- NormalizeData(obj, verbose = TRUE)
  obj <- FindVariableFeatures(obj, nfeatures = nfeatures, verbose = TRUE)

  obj <- ScaleData(
    obj,
    features   = VariableFeatures(obj),
    block.size = 200,
    verbose    = TRUE
  )

  obj <- RunPCA(obj,
                features = VariableFeatures(obj),
                npcs     = max(PCs_use),
                verbose  = TRUE)

  obj <- FindNeighbors(obj, dims = PCs_use)
  obj <- FindClusters(obj, resolution = 0.1)

  sweep.res.list <- paramSweep(obj, PCs = PCs_use, sct = FALSE)
  sweep.stats    <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn          <- find.pK(sweep.stats)

  pk_best <- bcmvn %>%
    dplyr::arrange(desc(BCmetric)) %>%
    dplyr::pull(pK) %>%
    .[1] %>%
    as.character() %>%
    as.numeric()

  message("样本 ", sample_id, " 最佳 pK = ", pk_best)

  annotations    <- obj$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  message("homotypic proportion = ", round(homotypic.prop, 3))

  nExp_poi      <- round(exp_rate * ncol(obj))
  nExp_poi.adj  <- round(nExp_poi * (1 - homotypic.prop))

  message("预期双细胞总数（未校正）=", nExp_poi,
          "；校正后 =", nExp_poi.adj)

  obj <- doubletFinder(
    obj,
    PCs       = PCs_use,
    pN        = pN,
    pK        = pk_best,
    nExp      = nExp_poi.adj,
    reuse.pANN = NULL,
    sct       = FALSE
  )

  meta_cols <- colnames(obj@meta.data)
  last_two  <- tail(meta_cols, 2)

  obj$Double_score <- obj@meta.data[[ last_two[1] ]]
  obj$Is_Double    <- obj@meta.data[[ last_two[2] ]]

  message("DoubletFinder 列名：", paste(last_two, collapse = ", "))
  message("Singlet / Doublet 计数：")
  print(table(obj$Is_Double))

  obj_singlet <- subset(obj, subset = Is_Double == "Singlet")

  message("过滤后细胞数：", ncol(obj_singlet))
  message("===== Done sample: ", sample_id, " =====")

  return(obj_singlet)
}

singlet_list <- lapply(seurat_list, run_doubletfinder_one)

rm(seurat_list)
gc()

seurat_singlet <- singlet_list[[1]]
if (length(singlet_list) > 1) {
  seurat_singlet <- merge(
    seurat_singlet,
    y       = singlet_list[-1],
    project = "OA2_nodouble"
  )
}

rm(singlet_list)
gc()

message("合并后 Singlet 总细胞数：", ncol(seurat_singlet))

seurat_singlet <- seurat_singlet %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 3000)

seurat_singlet <- seurat_singlet %>%
  ScaleData(
    features   = VariableFeatures(.),
    block.size = 200,
    verbose    = TRUE
  )

seurat_singlet <- seurat_singlet %>%
  RunPCA(features = VariableFeatures(.), npcs = 30) %>%
  RunUMAP(dims = 1:30) %>%
  RunTSNE(dims = 1:30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.1)

saveRDS(seurat_singlet, file = "nodouble.rds")

message(">>> 全部完成，结果已保存到 nodouble.rds")
