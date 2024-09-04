requireNamespace("jsonlite", quietly = TRUE)
requireNamespace("yaml", quietly = TRUE)
library(purrr, warn.conflicts = FALSE)
library(rlang, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input = "output/temp/task_info.yaml",
  output = "output/test/task_info.json"
)
## VIASH END

info <- yaml::yaml.load_file(par$input)
# ↑ this could be used as the new format

# construct v1 format
out <- list(
  task_id = info$name,
  commit_sha = NA_character_,
  task_name = info$label,
  task_summary = info$summary,
  task_description = paste0(info$motivation, "\n\n", info$description),
  repo = "openproblems-bio/openproblems-v2",
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
