### VIASH START
par <- list(
  input = "resources_test/openproblems/task_results_v4/raw/task_info.yaml",
  output = "task_info.json"
)
## VIASH END

source(file.path(meta$resources_dir, "functions.R"))

cat("====== Get task info ======\n")

`%||%` <- rlang::`%||%`
cat("\n>>> Reading input files...\n")
cat("Reading task info from '", par$input, "'...\n", sep = "")
task_info_yaml <- yaml::read_yaml(par$input)

cat("\n>>> Getting references...\n")
bibliography <- read_bibliography(
  file.path(meta$resources_dir, "bibliography.bib")
)
references <- get_references_list(task_info_yaml$references, bibliography)
str(references)

cat("\n>>> Getting authors...\n")
authors <- get_authors_list(task_info_yaml$authors)
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

cat("\n>>> Validating output against schema...\n")
results_schemas <- file.path(meta$resources_dir, "schemas", "results_v4")
ajv_args <- paste(
  "validate",
  "--spec draft2020",
  "-s",
  file.path(results_schemas, "task_info.json"),
  "-r",
  file.path(results_schemas, "core.json"),
  "-d",
  par$output
)

cat("Running validation command:", "ajv", ajv_args, "\n")
cat("Output:\n")
validation_result <- system2("ajv", ajv_args)

if (validation_result == 0) {
  cat("JSON validation passed successfully!\n")
} else {
  cat("JSON validation failed!\n")
  stop("Output JSON does not conform to schema")
}

cat("\n>>> Done!\n")
