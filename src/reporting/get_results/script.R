## VIASH START
raw_dir <- "resources_test/openproblems/task_results_v4/raw"
processed_dir <- "resources_test/openproblems/task_results_v4/processed"

par <- list(
  # Inputs
  input_scores = paste0(raw_dir, "/score_uns.yaml"),
  input_trace = paste0(raw_dir, "/trace.txt"),
  input_dataset_info = paste0(processed_dir, "/dataset_info.json"),
  input_method_info = paste0(processed_dir, "/method_info.json"),
  input_method_configs = paste0(raw_dir, "/method_configs.yaml"),
  input_metric_info = paste0(processed_dir, "/metric_info.json"),
  # Outputs
  output = paste0(processed_dir, "/results.json")
)
## VIASH END

################################################################################
#                               FUNCTIONS
################################################################################

parse_exit_code <- function(exit_codes) {
  as.integer(exit_codes)
}

parse_duration <- function(durations) {
  durations |>
    toupper() |>
    lubridate::duration() |>
    as.numeric()
}

parse_cpu_pct <- function(cpu_pcts) {
  cpu_pcts |>
    stringr::str_remove(" *%") |>
    as.numeric()
}

parse_memory <- function(memories) {
  values <- memories |>
    stringr::str_remove("[[:blank:][:alpha:]]+") |>
    as.numeric()

  units <- stringr::str_remove(memories, "[[:digit:]\\.[:blank:]]+")

  multipliers <- dplyr::case_when(
    units == "TB" ~ 1024 * 1024,
    units == "GB" ~ 1024,
    units == "MB" ~ 1,
    units == "KB" ~ 1 / 1024,
    units == "B" ~ 1 / 1024 / 1024,
    TRUE ~ NA
  )

  (values * multipliers) |>
    ceiling() |>
    as.integer()
}

missing_to_empty <- function(
  values,
  mode = c("character", "numeric", "integer")
) {
  mode <- match.arg(mode)

  if (is.null(values) || (length(values) == 1 && is.na(values))) {
    switch(
      mode,
      character = character(0),
      numeric = numeric(0),
      integer = integer(0)
    )
  } else {
    values
  }
}

map_missing_to_empty <- function(
  values_list,
  mode = c("character", "numeric")
) {
  purrr::map(values_list, missing_to_empty, mode = mode)
}

################################################################################
#                              MAIN SCRIPT
################################################################################

cat("====== Get results ======\n")

cat("\n>>> Reading input files...\n")
cat("Reading method info from '", par$input_method_info, "'...\n", sep = "")
method_info <- jsonlite::read_json(par$input_method_info)
cat("Reading dataset info from '", par$input_dataset_info, "'...\n", sep = "")
dataset_info <- jsonlite::read_json(par$input_dataset_info)
cat("Reading metric info from '", par$input_metric_info, "'...\n", sep = "")
metric_info <- jsonlite::read_json(par$input_metric_info)
cat("Reading scores from '", par$input_scores, "'...\n", sep = "")
scores <- yaml::yaml.load_file(par$input_scores) |>
  purrr::map_dfr(\(.x) {
    .x[c("dataset_id", "method_id", "metric_ids", "metric_values")] |>
      tibble::as_tibble()
  }) |>
  dplyr::rename(
    dataset_name = dataset_id,
    method_name = method_id,
    metric_name = metric_ids,
    metric_value = metric_values
  )
cat("Reading execution trace from '", par$input_trace, "'...\n", sep = "")
method_names <- purrr::map_chr(method_info, "name")
metric_components <- unique(purrr::map_chr(metric_info, "component_name"))
trace <- readr::read_tsv(
  par$input_trace,
  col_types = readr::cols(
    task_id = readr::col_integer(),
    submit = readr::col_datetime(),
    .default = readr::col_character(),
  ),
  na = c("", "-", "NA")
) |>
  # Only keep the most recent run of each process
  dplyr::group_by(name) |>
  dplyr::slice_max(submit) |>
  dplyr::ungroup() |>
  # Separate process name and id
  dplyr::mutate(name_copy = name) |>
  tidyr::separate_wider_delim(name_copy, " ", names = c("process", "id")) |>
  # Extract component from process name
  dplyr::mutate(
    component = purrr::map_chr(process, \(.process) {
      rev(stringr::str_split(.process, ":")[[1]])[1]
    })
  ) |>
  dplyr::mutate(component = stringr::str_remove(component, "_process")) |>
  # Only keep method and metric components
  dplyr::filter(
    component %in% method_names | component %in% metric_components
  ) |>
  dplyr::mutate(id = stringr::str_remove_all(id, "\\(|\\)")) |>
  # Split ID into dataset, method, metric
  tidyr::separate_wider_delim(
    id,
    delim = ".",
    names = c("dataset_name", "method_name", "metric_component"),
    too_few = "align_start"
  ) |>
  # Parse resources
  dplyr::mutate(
    run_exit_code = parse_exit_code(exit),
    run_duration_secs = parse_duration(realtime),
    run_cpu_pct = parse_cpu_pct(`%cpu`),
    run_peak_memory_mb = parse_memory(peak_vmem),
    run_disk_read_mb = parse_memory(rchar),
    run_disk_write_mb = parse_memory(wchar)
  ) |>
  # Select columns
  dplyr::select(
    name,
    process,
    component,
    dataset_name,
    method_name,
    metric_component,
    tidyselect::starts_with("run_")
  )

# Dataset names in the trace may have normalisations appended, map back to the name
dataset_names <- purrr::map_chr(dataset_info, "name")
process_datasets <- unique(trace$dataset_name)
dataset_map <- purrr::map_chr(process_datasets, function(.dataset) {
  dataset_names[stringr::str_detect(.dataset, dataset_names)][1]
}) |>
  purrr::set_names(process_datasets)
trace$dataset_name <- dataset_map[trace$dataset_name]

cat("\n>>> Extracting resources...\n")
cat("Extracting method resources...\n", sep = "")
method_resources <- trace |>
  dplyr::filter(component %in% method_names) |>
  dplyr::group_by(dataset_name, method_name) |>
  dplyr::summarise(
    run_exit_code = list(run_exit_code),
    run_duration_secs = list(run_duration_secs),
    run_cpu_pct = list(run_cpu_pct),
    run_peak_memory_mb = list(run_peak_memory_mb),
    run_disk_read_mb = list(run_disk_read_mb),
    run_disk_write_mb = list(run_disk_write_mb),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    succeeded = purrr::map_lgl(run_exit_code, ~ all(.x == 0)),
    run_exit_code = map_missing_to_empty(run_exit_code, mode = "integer"),
    run_duration_secs = map_missing_to_empty(
      run_duration_secs,
      mode = "numeric"
    ),
    run_cpu_pct = map_missing_to_empty(run_cpu_pct, mode = "numeric"),
    run_peak_memory_mb = map_missing_to_empty(
      run_peak_memory_mb,
      mode = "numeric"
    ),
    run_disk_read_mb = map_missing_to_empty(run_disk_read_mb, mode = "numeric"),
    run_disk_write_mb = map_missing_to_empty(
      run_disk_write_mb,
      mode = "numeric"
    )
  ) |>
  dplyr::relocate(succeeded, .after = method_name)

cat("Extracting metric resources...\n", sep = "")
metric_resources <- trace |>
  dplyr::filter(component %in% metric_components) |>
  dplyr::group_by(dataset_name, method_name, metric_component) |>
  dplyr::summarise(
    run_exit_code = list(run_exit_code),
    run_duration_secs = list(run_duration_secs),
    run_cpu_pct = list(run_cpu_pct),
    run_peak_memory_mb = list(run_peak_memory_mb),
    run_disk_read_mb = list(run_disk_read_mb),
    run_disk_write_mb = list(run_disk_write_mb),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    succeeded = purrr::map_lgl(run_exit_code, ~ all(.x == 0)),
    run_exit_code = map_missing_to_empty(run_exit_code, mode = "integer"),
    run_duration_secs = map_missing_to_empty(
      run_duration_secs,
      mode = "numeric"
    ),
    run_cpu_pct = map_missing_to_empty(run_cpu_pct, mode = "numeric"),
    run_peak_memory_mb = map_missing_to_empty(
      run_peak_memory_mb,
      mode = "numeric"
    ),
    run_disk_read_mb = map_missing_to_empty(run_disk_read_mb, mode = "numeric"),
    run_disk_write_mb = map_missing_to_empty(
      run_disk_write_mb,
      mode = "numeric"
    )
  ) |>
  dplyr::relocate(succeeded, .after = method_name)

cat("\n>>> Summarising results...\n")
metric_component_names <- purrr::map_chr(metric_info, "component_name")
metric_component_map <- purrr::map_chr(metric_info, "name") |>
  purrr::set_names(metric_component_names)
results <- scores |>
  # There shouldn't be any but remove missing/NaN values just in case
  dplyr::filter(
    !is.na(metric_value) & is.finite(metric_value)
  ) |>
  dplyr::arrange(dataset_name, method_name, metric_name) |>
  dplyr::group_by(dataset_name, method_name) |>
  dplyr::summarise(
    metric_names = list(metric_name),
    metric_values = list(metric_value),
    .groups = "drop"
  ) |>
  dplyr::full_join(method_resources, by = c("dataset_name", "method_name")) |>
  dplyr::mutate(
    metric_components = purrr::map2(
      dataset_name,
      method_name,
      function(.dataset, .method) {
        metric_resources |>
          dplyr::filter(
            dataset_name == .dataset,
            method_name == .method
          ) |>
          dplyr::mutate(
            metric_names = purrr::map(metric_component, function(.component) {
              metric_component_map[names(metric_component_map) == .component]
            })
          ) |>
          dplyr::select(
            component_name = metric_component,
            metric_names,
            succeeded,
            tidyselect::starts_with("run_")
          )
      }
    )
  ) |>
  # TODO: Add these once available in output
  dplyr::mutate(
    paramset_name = NA,
    paramset = NA
  ) |>
  dplyr::mutate(
    metric_names = map_missing_to_empty(metric_names, mode = "character"),
    metric_values = map_missing_to_empty(metric_values, mode = "numeric")
  ) |>
  dplyr::select(
    dataset_name,
    method_name,
    paramset_name,
    paramset,
    succeeded,
    tidyselect::starts_with("run_"),
    metric_names,
    metric_values,
    metric_components
  )

dplyr::glimpse(results)

cat("\n>>> Writing output files...\n")
cat("Writing results to '", par$output, "'...\n", sep = "")
jsonlite::write_json(
  results,
  par$output,
  pretty = TRUE,
  null = "null",
  na = "null"
)

cat("\n>>> Validating output against schema...\n")
ajv_args <- paste(
  "validate",
  "--spec draft2020",
  "-s", file.path(meta$resources_dir, "schemas", "results_schema.json"),
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
