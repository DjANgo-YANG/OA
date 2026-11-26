rm(list=ls())
gc()

source(file.path(Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."), "R", "project_paths.R"))

library(Seurat)
work_dir <- project_path("results_root", "single_cell/merged")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)
pathSamples=project_path("data_root", "single_cell/10x")

samples <- list.files(pathSamples)
seurat_list <- list()
for (sample in samples){

  data.path <- file.path(pathSamples, sample)

  seurat_data <- Read10X(data.dir = data.path)

  seurat_obj <- CreateSeuratObject(counts = seurat_data, project = sample,
                                   min.features = 200,
                                   min.cells = 0)

  seurat_list <- append(seurat_list, seurat_obj)
}

seurat_object <- merge(seurat_list[[1]], y = seurat_list[-1],
                       add.cell.ids = samples)

seurat_object$sample <- Idents(seurat_object)

seurat_object <- JoinLayers(seurat_object)
table(Idents(seurat_object))

print(seurat_object)
dim(seurat_object)

saveRDS(seurat_object,file = "1.rds")
