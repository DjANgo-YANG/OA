env_path <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

project_root <- env_path(
  "OAGLY_PROJECT_ROOT",
  normalizePath(".", winslash = "/", mustWork = FALSE)
)
data_root <- env_path("OAGLY_DATA_ROOT", file.path(project_root, "data"))
results_root <- env_path("OAGLY_RESULTS_ROOT", file.path(project_root, "results"))

paths <- list(
  project_root = project_root,
  data_root = data_root,
  raw_data = file.path(data_root, "raw"),
  metadata = file.path(data_root, "metadata"),
  references = env_path("OAGLY_REFERENCE_ROOT", file.path(data_root, "references")),
  derived_data = file.path(data_root, "derived"),
  results_root = results_root,
  figures = file.path(results_root, "figures"),
  tables = file.path(results_root, "tables"),
  objects = file.path(results_root, "objects")
)

invisible(lapply(paths[c("results_root", "figures", "tables", "objects")], dir.create,
  recursive = TRUE, showWarnings = FALSE
))
