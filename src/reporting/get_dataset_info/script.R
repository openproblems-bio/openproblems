## VIASH START
par <- list(
  input = "resources_test/openproblems/task_results_v4/raw/dataset_uns.yaml",
  output = "resources_test/openproblems/task_results_v4/processed/dataset_info.json"
)
## VIASH END

source(file.path(meta$resources_dir, "functions.R"))

`%||%` <- rlang::`%||%`

cat("====== Get dataset info ======\n")

cat("\n>>> Reading input files...\n")
cat("Reading dataset uns from '", par$input, "'...\n", sep = "")
dataset_uns <- yaml::yaml.load_file(
  par$input,
  # Read file sizes as floats to avoid issues with big integers
  handlers = list(int = \(x) {
    as.numeric(x)
  })
)

cat(
  "\n>>> Processing ",
  length(dataset_uns),
  " datasets...\n",
  sep = ""
)
bibliography <- read_bibliography(
  file.path(meta$resources_dir, "bibliography.bib")
)
dataset_info_json <- purrr::map(dataset_uns, function(.dataset) {
  cat("Processing dataset uns '", .dataset$dataset_id, "'\n", sep = "")

  authors <- get_authors_list(.dataset$authors)

  if ("dataset_reference" %in% names(.dataset)) {
    reference_name <- "dataset_reference"
  } else if ("data_reference" %in% names(.dataset)) {
    reference_name <- "data_reference"
  } else {
    warning(
      "No reference found in dataset uns for '", .dataset$dataset_id, "'",
      immediate. = TRUE
    )
    reference_name <- "NO_REFERENCE"
  }

  references <- get_references_list(.dataset[[reference_name]], bibliography)

  if ("dataset_url" %in% names(.dataset)) {
    url_name <- "dataset_url"
  } else if ("data_url" %in% names(.dataset)) {
    url_name <- "data_url"
  } else {
    warning(
      "No URL found in dataset uns for '", .dataset$dataset_id, "'",
      immediate. = TRUE
    )
    url_name <- "NO_URL"
  }

  list(
    name = jsonlite::unbox(.dataset$dataset_id),
    label = jsonlite::unbox(.dataset$dataset_name %||% .dataset$dataset_id),
    commit = jsonlite::unbox(.dataset$dataset_commit %||% "missing-sha"),
    summary = .dataset$dataset_summary %||%
        .dataset$dataset_description %||%
        .dataset$dataset_id |>
      stringr::str_trim() |>
      stringr::str_remove_all('(^"|"$|^\'|\'$)') |>
      jsonlite::unbox(),
    description = .dataset$dataset_description %||%
        .dataset$dataset_summary %||%
        .dataset$dataset_id |>
      stringr::str_trim() |>
      stringr::str_remove_all('(^"|"$|^\'|\'$)') |>
      jsonlite::unbox(),
    source_urls = .dataset[[url_name]] %||% character(0),
    common_dataset_names = .dataset$common_dataset_id,
    modalities = .dataset$dataset_modality %||%
      .dataset$modality %||%
      character(0),
    organisms = .dataset$dataset_organism %||% character(0),
    authors = authors,
    references = references,
    date_created = jsonlite::unbox(.dataset$date_created),
    file_size_mb = if (is.null(.dataset$file_size)) {
      NULL
    } else {
      jsonlite::unbox(.dataset$file_size / 1048576)
    }
  )
})

dataset_names <- purrr::map_chr(dataset_info_json, "name")
if (any(duplicated(dataset_names))) {
  warning(
    "Duplicate dataset names found: ",
    paste(dataset_names[duplicated(dataset_names)], collapse = ", "),
    "\nOnly the first instance will be kept",
    immediate. = TRUE,
    call. = FALSE
  )

  dataset_info_json <- dataset_info_json[!duplicated(dataset_names)]
}

cat("\n>>> Writing output files...\n")
cat("Writing dataset info to '", par$output, "'...\n", sep = "")
jsonlite::write_json(
  dataset_info_json,
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
  file.path(results_schemas, "dataset_info.json"),
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
