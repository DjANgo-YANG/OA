project_root <- normalizePath(
  Sys.getenv("OAGLY_PROJECT_ROOT", unset = "."),
  winslash = "/",
  mustWork = FALSE
)

local_config <- file.path(project_root, "config", "paths.R")
example_config <- file.path(project_root, "config", "paths.example.R")
config_file <- Sys.getenv(
  "OAGLY_PATH_CONFIG",
  unset = if (file.exists(local_config)) local_config else example_config
)

if (!file.exists(config_file)) {
  stop("Path configuration was not found: ", config_file)
}

source(config_file, local = .GlobalEnv, encoding = "UTF-8")

path_overrides <- setNames(character(), character())
override_file <- file.path(project_root, "config", "path_overrides.local.R")
if (file.exists(override_file) && normalizePath(config_file, winslash = "/") ==
    normalizePath(local_config, winslash = "/", mustWork = FALSE)) {
  source(override_file, local = .GlobalEnv, encoding = "UTF-8")
}

project_path <- function(key, ...) {
  base <- paths[[key]]
  if (is.null(base) || !nzchar(base)) {
    stop("Unknown or empty path key: ", key)
  }
  override <- unname(path_overrides[file.path(key, ...)])
  if (!is.na(override)) return(override)
  file.path(base, ...)
}
