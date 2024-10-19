requireNamespace("jsonlite", quietly = TRUE)
requireNamespace("yaml", quietly = TRUE)
library(purrr, warn.conflicts = FALSE)
library(rlang, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input = "resources_test/openproblems/task_results_v3/raw/task_info.yaml",
  output = "output/test/task_info.json"
)
## VIASH END

info <- yaml::yaml.load_file(par$input)
# ↑ this could be used as the new format

# construct v1 format
repo <-
  if ("name" %in% names(info) && "organization" %in% names(info)) {
    paste0(info$organization, "/", info$name)
  } else {
    "openproblems-bio/openproblems"
  }
description <-
  if ("motivation" %in% names(info)) {
    paste0(info$motivation, "\n\n", info$description)
  } else {
    info$description
  }
out <- list(
  task_id = info$name,
  commit_sha = NA_character_,
  task_name = info$label,
  task_summary = info$summary,
  task_description = description,
  repo = repo,
  authors = info$authors
)

# show warning when certain data is missing and return null?
for (n in names(out)) {
  if (is.null(out[[n]])) {
    out_as_str <- jsonlite::toJSON(out, auto_unbox = TRUE, pretty = TRUE)
    stop("missing value for value '", n, "' in ", out_as_str)
  }
}

jsonlite::write_json(
  out,
  par$output,
  auto_unbox = TRUE,
  pretty = TRUE
)
