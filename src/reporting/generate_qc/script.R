## VIASH START
processed_dir <- "resources_test/openproblems/task_results_v4/processed"

par <- list(
  # Inputs
  input_task_info = paste0(processed_dir, "/task_info.json"),
  input_method_info = paste0(processed_dir, "/method_info.json"),
  input_metric_info = paste0(processed_dir, "/metric_info.json"),
  input_dataset_info = paste0(processed_dir, "/dataset_info.json"),
  input_results = paste0(processed_dir, "/results.json"),
  # Outputs
  output = "quality_control.json"
)
## VIASH END

################################################################################
#                               FUNCTIONS
################################################################################

create_qc_entry <- function(
  category,
  label,
  value,
  severity_value,
  condition,
  message
) {
  # If values are missing, set to -1
  # This can happen if a method/metric is not run and therefore has no results
  if (is.null(value) || is.na(value) || length(value) == 0) {
    value <- -1
  }

  if (is.null(severity_value) || is.na(severity_value) || length(severity_value) == 0) {
    severity_value <- -1
  }

  severity <- dplyr::case_when(
    severity_value < 0 ~ 3L,
    severity_value < 1 ~ 0L,
    severity_value < 2 ~ 1L,
    severity_value < 3 ~ 2L,
    TRUE ~ 3L
  )

  list(
    category = category,
    label = label,
    value = value,
    severity = severity,
    severity_value = severity_value,
    condition = condition,
    message = message
  )
}

percent_missing <- function(items, field) {
  is_missing <- purrr::map_lgl(items, \(.item) {
    if (field == "references") {
      return(references_missing(.item))
    }

    field_value <- .item[[field]]
    is.null(field_value) ||
      is.na(field_value) ||
      (is.character(field_value) && field_value == "")
  })

  mean(is_missing)
}

references_missing <- function(item) {
  # Special case for control methods without references
  if ("type" %in% names(item) && item$type == "control_method") {
    return(FALSE)
  }

  references <- item$references
  if (length(references) == 0) {
    return(TRUE)
  }

  if (
    length(references) == 2 && all(c("doi", "bibtex") %in% names(references))
  ) {
    if (length(references$doi) == 0 && length(references$bibtex) == 0) {
      return(TRUE)
    }
  }

  return(FALSE)
}

check_info_fields <- function(info, type, expected_fields, task_name) {
  category <- paste(stringr::str_to_title(type), "info")

  purrr::map(expected_fields, function(.field) {
    pct_missing <- percent_missing(info, .field)
    create_qc_entry(
      category = category,
      label = paste0("Info field '", .field, "' % missing"),
      value = pct_missing,
      severity_value = ifelse(pct_missing > 0, 3.0, 0.0),
      condition = "pct_missing <= 0",
      message = paste0(
        category, " field '", .field, "' should be defined\n",
        "  Task: ", task_name, "\n",
        "  Field: ", .field, "\n",
        "  Percentage missing: ", round(pct_missing * 100, 0)
      )
    )
  })
}

check_missing_results <- function(
  results_long,
  name,
  type,
  n_datasets,
  n_methods,
  n_metrics,
  task_name
) {
  n_expected <- switch(
    type,
    "dataset" = n_methods * n_metrics,
    "method" = n_datasets * n_metrics,
    "metric" = n_datasets * n_methods
  )

  name_col <- paste0(type, "_name")
  n_results <- results_long |>
    dplyr::filter(!!rlang::sym(name_col) == name) |>
    nrow()
  pct_missing <- 1 - (n_results / n_expected)

  title <- type |>
    stringr::str_replace_all("_", " ") |>
    stringr::str_to_sentence()

  create_qc_entry(
    category = "Raw results",
    label = paste0(title, " '", name, "' % missing"),
    value = pct_missing,
    severity_value = pct_missing / 0.1,
    condition = "pct_missing <= 0.1",
    message = paste0(
      "Percentage of missing results should be less than 10%\n",
      "  Task: ", task_name, "\n",
      "  ", title, ": ", name, "\n",
      "  Number of results: ", n_results, "\n",
      "  Expected number of results: ", n_expected, "\n",
      "  Percentage missing: ", round(pct_missing * 100, 0), "%\n"
    )
  )
}

check_failed_processes <- function(results, name, type, task_name) {
  name_col <- paste0(type, "_name")
  results_name <- results |>
    dplyr::filter(!!rlang::sym(name_col) == name)

  n_expected <- nrow(results_name)
  n_succeeded <- sum(results_name$succeeded)
  pct_failed <- 1 - (n_succeeded / n_expected)

  title <- type |>
    stringr::str_replace_all("_", " ") |>
    stringr::str_to_sentence()

  create_qc_entry(
    category = "Raw results",
    label = paste0(title, " '", name, "' % failed"),
    value = pct_failed,
    severity_value = pct_failed / 0.1,
    condition = "pct_failed <= 0.1",
    message = paste0(
      "Percentage of failed processes should be less than 10%\n",
      "  Task: ", task_name, "\n",
      "  ", title, ": ", name, "\n",
      "  Succeeded processes: ", n_succeeded, "\n",
      "  Attempted processes: ", n_expected, "\n",
      "  Percentage failed: ", round(pct_failed * 100, 0), "%\n"
    )
  )
}

check_metric_scaling <- function(
  results_long,
  metric,
  control_methods,
  task_name
) {
  `%||%` <- rlang::`%||%`

  metric_results <- results_long |>
    dplyr::filter(metric_name == metric) |>
    dplyr::select(-metric_name)

  if (
    nrow(metric_results) == 0 ||
      !any(control_methods %in% metric_results$method_name)
  ) {
    return(list())
  }

  control_range <- metric_results |>
    dplyr::filter(
      method_name %in% control_methods
    ) |>
    dplyr::group_by(dataset_name) |>
    dplyr::summarise(
      control_min = min(metric_value),
      control_max = max(metric_value)
    )

  scaled_metrics <- metric_results |>
    dplyr::left_join(control_range, by = "dataset_name") |>
    dplyr::mutate(
      scaled_value = (metric_value - control_min) / (control_max - control_min),
      outside = scaled_value < 0 | scaled_value > 1,
      pct_outside = dplyr::case_when(
        scaled_value < 0 ~ 0 - scaled_value,
        scaled_value > 1 ~ scaled_value - 1,
        TRUE ~ NA
      )
    )

  pct_outside <- sum(scaled_metrics$outside) / nrow(scaled_metrics)
  worst_score <- min(scaled_metrics$scaled_value)
  worst_pct_outside <- if (worst_score < 0) {
    max(scaled_metrics$pct_outside[scaled_metrics$scaled_value < 0])
  } else {
    0
  }
  best_score <- max(scaled_metrics$scaled_value)
  best_pct_outside <- if (best_score > 1) {
    max(scaled_metrics$pct_outside[scaled_metrics$scaled_value > 1])
  } else {
    0
  }

  metric_checks <- list(
    create_qc_entry(
      category = "Scaling",
      label = paste0("Metric '", metric, "' % outside range"),
      value = pct_outside,
      severity_value = pct_outside / 0.1,
      condition = "pct_outside <= 0.1",
      message = paste0(
        "Percentage of scaled scores outside control range should be less than 10%\n",
        "  Task: ", task_name, "\n",
        "  Metric: ", metric, "\n",
        "  Inside range: ", sum(!scaled_metrics$outside), "\n",
        "  Scaled scores: ", nrow(scaled_metrics), "\n",
        "  Percentage outside: ", round(pct_outside * 100, 0), "%\n"
      )
    ),
    create_qc_entry(
      category = "Scaling",
      label = paste0("Metric '", metric, "' worst score % outside range"),
      value = worst_pct_outside,
      severity_value = worst_pct_outside / 0.1,
      condition = "worst_pct_outside <= 0.1",
      message = paste0(
        "The worst scaled score should be less than 10% outside the control range\n",
        "  Task: ", task_name, "\n",
        "  Metric: ", metric, "\n",
        "  Worst score: ", worst_score, "\n",
        "  Percentage outside range: ", round(worst_pct_outside * 100, 0), "%\n"
      )
    ),
    create_qc_entry(
      category = "Scaling",
      label = paste0("Metric '", metric, "' best score % outside range"),
      value = best_pct_outside,
      severity_value = best_pct_outside / 0.1,
      condition = "best_pct_outside <= 0.1",
      message = paste0(
        "The best scaled score should be less than 10% outside the control range\n",
        "  Task: ", task_name, "\n",
        "  Metric: ", metric, "\n",
        "  Best score: ", best_score, "\n",
        "  Percentage outside range: ", round(best_pct_outside * 100, 0), "%\n"
      )
    )
  )

  method_metric_checks <- purrr::map(
    sort(unique(scaled_metrics$method_name)),
    function(.method) {
      check_method_metric_scaling(scaled_metrics, .method, task_name, metric)
    }
  ) |>
    purrr::list_flatten()

  c(metric_checks, method_metric_checks)
}

check_method_metric_scaling <- function(
  scaled_metrics,
  method,
  task_name,
  metric_name
) {
  method_scaled_metrics <- scaled_metrics |>
    dplyr::filter(method_name == method)

  worst_score <- min(method_scaled_metrics$scaled_value)
  worst_pct_outside <- if (worst_score < 0) {
    max(method_scaled_metrics$pct_outside[
      method_scaled_metrics$scaled_value < 0
    ])
  } else {
    0
  }
  best_score <- max(method_scaled_metrics$scaled_value)
  best_pct_outside <- if (best_score > 1) {
    max(method_scaled_metrics$pct_outside[
      method_scaled_metrics$scaled_value > 1
    ])
  } else {
    0
  }

  list(
    create_qc_entry(
      category = "Scaling",
      label = paste0("Worst '", metric_name, "' score for '", method, "'"),
      value = worst_score,
      severity_value = ifelse(worst_score < -1, worst_pct_outside, 0),
      condition = "worst_score < -1",
      message = paste0(
        "Method '", method,"' performs much worse than controls for metric' ", metric_name, "'\n",
        "  Task: ",task_name, "\n",
        "  Method: ", method, "\n",
        "  Metric: ", metric_name, "\n",
        "  Worst score: ", worst_score, "\n",
        "  Percentage outside range: ", round(worst_pct_outside * 100, 0), "%\n"
      )
    ),
    create_qc_entry(
      category = "Scaling",
      label = paste0("Best '", metric_name, "' score for '", method, "'"),
      value = best_score,
      severity_value = ifelse(best_score > 2, best_pct_outside, 0),
      condition = "best_score > 2",
      message = paste0(
        "Method '", method, "' performs much better than controls for metric '", metric_name, "'\n",
        "  Task: ", task_name, "\n",
        "  Method: ", method, "\n",
        "  Metric: ", metric_name, "\n",
        "  Best score: ", best_score, "\n",
        "  Percentage outside range: ", round(best_pct_outside * 100, 0), "%\n"
      )
    )
  )
}

################################################################################
#                              MAIN SCRIPT
################################################################################

cat("====== Generate QC ======\n")

cat("\n>>> Reading input files...\n")
cat("Reading task info from '", par$input_task_info, "'...\n", sep = "")
task_info <- jsonlite::read_json(par$input_task_info)

cat("Reading dataset info from '", par$input_dataset_info, "'...\n", sep = "")
dataset_info <- jsonlite::read_json(par$input_dataset_info)

cat("Reading method info from '", par$input_method_info, "'...\n", sep = "")
method_info <- jsonlite::read_json(par$input_method_info)

cat("Reading metric info from '", par$input_metric_info, "'...\n", sep = "")
metric_info <- jsonlite::read_json(par$input_metric_info)

cat("Reading results from '", par$input_results, "'...\n", sep = "")
results <- jsonlite::read_json(par$input_results, simplifyVector = TRUE)

cat("\n>>> Checking expected info fields...\n")

expected_task_fields <- c("name", "label", "summary", "description")
expected_dataset_fields <- c(
  "name",
  "label",
  "summary",
  "description",
  "references"
)
expected_method_fields <- c(
  "name",
  "label",
  "commit",
  "summary",
  "description",
  "references"
)
expected_metric_fields <- c(
  "name",
  "label",
  "commit",
  "summary",
  "description",
  "references"
)

task_name <- task_info$name %||% "unknown"

info_task <- check_info_fields(
  list(task_info),
  "task",
  expected_task_fields,
  task_name
)
info_datasets <- check_info_fields(
  dataset_info,
  "dataset",
  expected_dataset_fields,
  task_name
)
info_methods <- check_info_fields(
  method_info,
  "method",
  expected_method_fields,
  task_name
)
info_metrics <- check_info_fields(
  metric_info,
  "metric",
  expected_metric_fields,
  task_name
)

cat("\n>>> Checking missing results...\n")
results_long <- results |>
  dplyr::select(dataset_name, method_name, metric_names, metric_values) |>
  tidyr::unnest_longer(c("metric_names", "metric_values")) |>
  dplyr::rename(
    metric_name = metric_names,
    metric_value = metric_values
  ) |>
  dplyr::filter(!is.na(metric_value))

dataset_names <- purrr::map_chr(dataset_info, "name")
method_names <- purrr::map_chr(method_info, "name")
metric_names <- purrr::map_chr(metric_info, "name")

n_datasets <- length(dataset_names)
n_methods <- length(method_names)
n_metrics <- length(metric_names)

n_results_expected <- n_datasets * n_methods * n_metrics
n_results <- nrow(results_long)
pct_results_missing <- 1 - (n_results / n_results_expected)

results_task <- list(
  create_qc_entry(
    category = "Raw results",
    label = "Task number of results",
    value = n_results,
    severity_value = pct_results_missing / 0.1,
    condition = "length(dataset_info) * length(results) == length(method_info) * length(metric_info)",
    message = paste0(
      "Number of results should be equal to #datasets × #methods × #metrics \n",
      "  Task: ", task_name, "\n",
      "  Number of results: ", n_results, "\n",
      "  Number of datasets: ", n_datasets, "\n",
      "  Number of methods: ", n_methods, "\n",
      "  Number of metrics: ", n_metrics, "\n",
      "  Expected number of results: ", n_results_expected, "\n"
    )
  )
)

results_datasets <- purrr::map(dataset_names, function(.dataset) {
  check_missing_results(
    results_long,
    .dataset,
    "dataset",
    n_datasets,
    n_methods,
    n_metrics,
    task_name
  )
})
results_methods <- purrr::map(method_names, function(.method) {
  check_missing_results(
    results_long,
    .method,
    "method",
    n_datasets,
    n_methods,
    n_metrics,
    task_name
  )
})
results_metrics <- purrr::map(metric_names, function(.metric) {
  check_missing_results(
    results_long,
    .metric,
    "metric",
    n_datasets,
    n_methods,
    n_metrics,
    task_name
  )
})

cat("\n>>> Checking failed processes\n")
metric_component_results <- results |>
  dplyr::select(dataset_name, method_name, metric_components) |>
  tidyr::unnest(metric_components) |>
  dplyr::rename(metric_component_name = component_name)

n_processes <- nrow(results) + nrow(metric_component_results)
n_succeeded <- sum(results$succeeded) + sum(metric_component_results$succeeded)
pct_failed <- 1 - (n_succeeded / n_processes)

failed_task <- list(
  create_qc_entry(
    category = "Raw results",
    label = "Task number of successful processes",
    value = n_succeeded,
    severity_value = pct_failed / 0.1,
    condition = "sum(results$succeeded) + sum(metric_component_results$succeeded) == nrow(results) + nrow(metric_component_results)",
    message = paste0(
      "Number of successful processes should be equal to the number of attempted processes\n",
      "  Task: ", task_name, "\n",
      "  Succeeded processes: ", n_succeeded, "\n",
      "  Attempted processes: ", n_processes, "\n",
      "  Percentage failed: ", round(pct_failed * 100, 0), "%\n"
    )
  )
)

failed_datasets <- purrr::map(dataset_names, function(.dataset) {
  check_failed_processes(results, .dataset, "dataset", task_name)
})
failed_methods <- purrr::map(method_names, function(.method) {
  check_failed_processes(results, .method, "method", task_name)
})

failed_metrics <- purrr::map(
  unique(metric_component_results$metric_component_name),
  function(.component) {
    check_failed_processes(
      metric_component_results,
      .component,
      "metric_component",
      task_name
    )
  }
)

cat("\n>>> Checking control methods...\n")
is_control <- purrr::map_lgl(method_info, \(.method) {
  .method$type == "control_method"
})
control_methods <- method_names[is_control]

dataset_controls <- results_long |>
  dplyr::filter(method_name %in% control_methods) |>
  dplyr::select(dataset_name, method_name) |>
  dplyr::distinct() |>
  dplyr::group_by(dataset_name) |>
  dplyr::count(name = "n_controls") |>
  dplyr::ungroup() |>
  dplyr::mutate(dataset_name = factor(dataset_name, levels = dataset_names)) |>
  tidyr::complete(dataset_name, fill = list(n_controls = 0))

controls_datasets <- purrr::map(seq_len(nrow(dataset_controls)), function(.idx) {
  dataset_name <- dataset_controls$dataset_name[.idx]
  n_controls <- dataset_controls$n_controls[.idx]

  create_qc_entry(
    category = "Raw results",
    label = "Dataset number of control methods",
    value = n_controls,
    severity_value = ifelse(n_controls != length(control_methods), 3, 0),
    condition = "n_controls != length(control_methods)",
    message = paste0(
      "Number of successful control methods for a dataset should equal the number of controls\n",
      "  Task: ", task_name, "\n",
      "  Succeeded control_methods: ", n_controls, "\n",
      "  Total control methods: ", length(control_methods), "\n",
      "  Percentage succeeded: ", round(n_controls / length(control_methods) * 100, 0), "%\n"
    )
  )
})

metric_controls <- results_long |>
  dplyr::filter(method_name %in% control_methods) |>
  dplyr::select(method_name, metric_name) |>
  dplyr::group_by(metric_name) |>
  dplyr::count(name = "n_controls") |>
  dplyr::ungroup() |>
  dplyr::mutate(metric_name = factor(metric_name, levels = metric_names)) |>
  tidyr::complete(metric_name, fill = list(n_controls = 0))

n_expected <- length(dataset_names) * length(control_methods)
controls_metrics <- purrr::map(seq_len(nrow(metric_controls)), function(.idx) {
  metric_name <- metric_controls$metric_name[.idx]
  n_controls <- metric_controls$n_controls[.idx]

  create_qc_entry(
    category = "Raw results",
    label = "Metric number of control methods",
    value = n_controls,
    severity_value = ifelse(n_controls != n_expected, 3, 0),
    condition = "n_controls != length(datasets) * length(control_methods)",
    message = paste0(
      "Number of metric scores for control methods should be equal to #datasets × #control_methods\n",
      "  Task: ", task_name, "\n",
      "  Control method scores: ", n_controls, "\n",
      "  Expected control method scores: ", n_expected, "\n",
      "  Percentage succeeded: ", round(n_controls / n_expected * 100, 0), "%\n"
    )
  )
})

cat("\n>>> Checking metric scaling...\n")
scaling <- purrr::map(metric_names, function(.metric) {
  check_metric_scaling(results_long, .metric, control_methods, task_name)
}) |>
  purrr::list_flatten()

cat("\n>>> Collecting QC results...\n")
qc_results <- c(
  info_task,
  info_datasets,
  info_methods,
  info_metrics,
  results_task,
  results_datasets,
  results_methods,
  results_metrics,
  failed_task,
  failed_datasets,
  failed_methods,
  failed_metrics,
  controls_datasets,
  controls_metrics,
  scaling
)

cat("\n>>> Writing output file...\n")
cat("Writing quality control to '", par$output, "'...\n", sep = "")
jsonlite::write_json(
  qc_results,
  par$output,
  pretty = TRUE,
  null = "null",
  na = "null",
  auto_unbox = TRUE
)

cat("\n>>> Done!\n")
