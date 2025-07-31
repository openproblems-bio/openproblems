## VIASH START
par <- list(
  input = "resources_test/openproblems/task_results_v4/raw/dataset_uns.yaml",
  output = "resources_test/openproblems/task_results_v4/processed/dataset_info.json"
)
## VIASH END

cat("====== Get dataset info ======\n")

`%||%` <- rlang::`%||%`

cat("\n>>> Reading input files...\n")
cat("Reading dataset uns from '", par$input, "'...\n", sep = "")
dataset_uns <- yaml::yaml.load_file(
  par$input,
  # Read file sizes as floats to avoid issues with big integers
  handlers = list(int = \(x) {as.numeric(x)})
)

cat(
  "\n>>> Processing ", length(dataset_uns), " datasets...\n",
  sep = ""
)
dataset_info_json <- purrr::map(dataset_uns, function(.dataset) {
  cat("Processing dataset uns '", .dataset$dataset_id, "'\n", sep = "")

  if ("dataset_reference" %in% names(.dataset)) {
    reference_name <- "dataset_reference"
  } else if ("data_reference" %in% names(.dataset)) {
    reference_name <- "data_reference"
  } else {
    stop("No reference found in dataset uns for '", .dataset$dataset_id, "'")
  }
  references <- if (is.list(.dataset[[reference_name]])) {
    list(
      doi = .dataset[[reference_name]]$doi %||% character(0),
      bibtex = .dataset[[reference_name]]$bibtex %||% character(0)
    )
  } else {
    reference <- .dataset[[reference_name]]

    if (startsWith(reference, "@")) {
      list(
        doi = character(0),
        bibtex = reference
      )
    } else if (startsWith(reference, "1")) {
      list(
        doi = reference,
        bibtex = character(0)
      )
    } else {
      reference
    }
  }

  if ("dataset_url" %in% names(.dataset)) {
    url_name <- "dataset_url"
  } else if ("data_url" %in% names(.dataset)) {
    url_name <- "data_url"
  } else {
    stop("No URL found in dataset uns for '", .dataset$dataset_id, "'")
  }

  list(
    name = jsonlite::unbox(.dataset$dataset_id),
    label = jsonlite::unbox(.dataset$dataset_name),
    commit = jsonlite::unbox(.dataset$dataset_commit %||% "missing-sha"),
    summary = .dataset$dataset_summary |>
      stringr::str_trim() |>
      stringr::str_remove_all('(^"|"$|^\'|\'$)') |>
      jsonlite::unbox(),
    description = .dataset$dataset_description |>
      stringr::str_trim() |>
      stringr::str_remove_all('(^"|"$|^\'|\'$)') |>
      jsonlite::unbox(),
    source_url = jsonlite::unbox(.dataset[[url_name]]),
    common_dataset_names = .dataset$common_dataset_id,
    modalities = jsonlite::unbox(.dataset$dataset_modality),
    organisms = .dataset$dataset_organism,
    references = references,
    date_created = jsonlite::unbox(.dataset$date_created),
    file_size_mb = jsonlite::unbox(.dataset$file_size / 1048576)
  )
})

cat("\n>>> Writing output files...\n")
cat("Writing dataset info to '", par$output, "'...\n", sep = "")
jsonlite::write_json(
  dataset_info_json,
  par$output,
  pretty = TRUE,
  null = "null"
)

cat("\n>>> Validating output against schema...\n")
ajv_args <- paste(
  "validate",
  "--spec draft2020",
  "-s", file.path(meta$resources_dir, "schemas", "dataset_info_schema.json"),
  "-r", file.path(meta$resources_dir, "schemas", "references_schema.json"),
  "-d", par$output
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
