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

  references <- if (is.list(.dataset$dataset_reference)) {
    list(
      doi = .dataset$dataset_reference$doi %||% character(0),
      bibtex = .dataset$dataset_reference$bibtex %||% character(0)
    )
  } else {
    .dataset$dataset_reference
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
    source_url = jsonlite::unbox(.dataset$dataset_url),
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

cat("\n>>> Done!\n")
