rm(list=ls())

library(requirements)

req_dir <- function(path = ".") {
  files <- dir(path, recursive = TRUE, include.dirs = FALSE, full.names = TRUE)
  flat_map_chr <- function(.x, .f, ...) {
      purrr::flatten_chr(purrr::map(.x, .f, ...))
    }
  sort(unique(flat_map_chr(files, req_file)))
}

req_dir(path = ".")

