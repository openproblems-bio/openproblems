### VIASH START
par <- list(
  input = "resources_test/openproblems/task_results_v4/raw/task_info.yaml",
  output = "task_info.json"
)
### VIASH END

cat("====== Get task info ======\n")

`%||%` <- rlang::`%||%`

cat("\n>>> Reading input files...\n")
cat("Reading task info from '", par$input, "'...\n", sep = "")
task_info_yaml <- yaml::read_yaml(par$input)

cat("\n>>> Getting references...\n")
references <- if (!is.null(task_info_yaml$references)) {
  list(
    doi = task_info_yaml$references$doi %||% character(0),
    bibtex = task_info_yaml$references$bibtex %||% character(0)
  )
} else {
  list(doi = character(0), bibtex = character(0))
}
str(references)

cat("\n>>> Getting authors...\n")
authors <- purrr::map(task_info_yaml$authors, function(.author) {
  other_fields <- setdiff(names(.author$info), c("github", "orcid"))

  list(
    name = jsonlite::unbox(.author$name),
    roles = .author$roles %||% character(0),
    github = jsonlite::unbox(.author$info$github),
    orcid = jsonlite::unbox(.author$info$orcid),
    info = .author$info[other_fields]
  )
})
cat("Found", length(authors), "authors\n")

cat("\n>>> Creating JSON list...\n")
task_info_json <- list(
  name = jsonlite::unbox(sub("^task_", "", task_info_yaml$name)), # Remove "task_" prefix
  commit = jsonlite::unbox(NA_character_), # TODO: Add when available in task_info.yaml
  label = jsonlite::unbox(task_info_yaml$label),
  summary = task_info_yaml$summary |>
    stringr::str_trim() |>
    stringr::str_remove_all('(^"|"$|^\'|\'$)') |>
    jsonlite::unbox(),
  description = task_info_yaml$description |>
    stringr::str_trim() |>
    stringr::str_remove_all('(^"|"$|^\'|\'$)') |>
    jsonlite::unbox(),
  repository = jsonlite::unbox(task_info_yaml$links$repository),
  authors = authors,
  license = jsonlite::unbox(task_info_yaml$license),
  references = references,
  version = jsonlite::unbox(task_info_yaml$version),
  is_prerelease = jsonlite::unbox(TRUE)
)
str(task_info_json)

cat("\n>>> Writing output files...\n")
cat("Writing task info to '", par$output, "'...\n", sep = "")
jsonlite::write_json(
  task_info_json,
  par$output,
  pretty = TRUE,
  null = "null"
)

cat("\n>>> Done!\n")
