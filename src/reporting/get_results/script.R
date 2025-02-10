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
    output_results = paste0(processed_dir, "/results.json"),
    output_metric_execution_info = paste0(
        processed_dir, "/metric_execution_info.json"
    )
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

################################################################################
#                              MAIN SCRIPT
################################################################################

cat("====== Get results ======\n")

cat("\n>>> Reading input files...\n")
cat("Reading method info from '", par$input_method_info, "'...\n", sep = "")
method_info <- jsonlite::read_json(par$input_method_info)
cat("Reading dataset info from '", par$input_dataset_info, "'...\n", sep = "")
dataset_info <- jsonlite::read_json(par$input_dataset_info)
cat("Reading scores from '", par$input_scores, "'...\n", sep = "")
scores <- yaml::yaml.load_file(par$input_scores) |>
    purrr::map_dfr(\(.x) {
        .x[c("dataset_id", "method_id", "metric_ids", "metric_values")] |>
            as.data.frame()
    }) |>
    dplyr::rename(
        dataset_name = dataset_id,
        method_name = method_id,
        metric_name = metric_ids,
        metric_value = metric_values
    )
cat("Reading execution trace from '", par$input_trace, "'...\n", sep = "")
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
    # Parse resources
    dplyr::mutate(
        run_exit_code = parse_exit_code(exit),
        run_duration_secs = parse_duration(realtime),
        run_cpu_pct = parse_cpu_pct(`%cpu`),
        run_peak_memory_mb = parse_memory(peak_vmem),
        run_disk_read_mb = parse_memory(rchar),
        run_disk_write_mb = parse_memory(wchar)
    )

cat("\n>>> Mapping names to processes...\n")
# Get names from info files instead of scores to avoid missing anything that
# didn't result in any metric scores (i.e. a method that always fails)
method_names <- purrr::map_chr(method_info, "name")
dataset_names <- purrr::map_chr(dataset_info, "name")

cat("Mapping methods...\n", sep = "")
method_processes <- purrr::map_dfr(method_names, \(.method) {
    is_method <- stringr::str_detect(
        trace$name,
        paste0(":", .method, ":")
    )

    processes <- unique(trace$name[is_method])

    data.frame(method_name = .method, process_name = processes)
})

cat("Mapping datasets...\n", sep = "")
dataset_processes <- purrr::map_dfr(dataset_names, \(.dataset) {
    is_dataset <- stringr::str_detect(
        trace$name,
        paste0("[(.]", .dataset, "[)./]")
    )

    processes <- unique(trace$name[is_dataset])

    data.frame(dataset_name = .dataset, process_name = processes)
})

cat("\n>>> Extracting resources...\n")
resources <- trace |>
    # Add method names
    dplyr::left_join(method_processes, by = c(name = "process_name")) |>
    # Add dataset names
    dplyr::left_join(dataset_processes, by = c(name = "process_name")) |>
    # Select only processes with a method name
    dplyr::filter(!is.na(method_name)) |>
    # Summarise the resources columns
    dplyr::group_by(dataset_name, method_name) |>
    dplyr::summarise(
        run_exit_code = list(run_exit_code),
        run_duration_secs = list(run_duration_secs),
        run_cpu_pct = list(run_cpu_pct),
        run_peak_memory_mb = list(run_peak_memory_mb),
        run_disk_read_mb = list(run_disk_read_mb),
        run_disk_write_mb = list(run_disk_write_mb),
        .groups = "drop"
    )

cat("\n>>> Summarising results...\n")
results <- scores |>
    # There shouldn't be any but remove missing/NaN values just in case
    dplyr::filter(
        !is.na(metric_value) & is.finite(metric_value)
    ) |>
    dplyr::group_by(dataset_name, method_name) |>
    dplyr::summarise(
        metric_names = list(metric_name),
        metric_values = list(metric_value),
        .groups = "drop"
    ) |>
    dplyr::left_join(resources, by = c("dataset_name", "method_name")) |>
    # TODO: Add these once available in output
    dplyr::mutate(
        paramset_name = NA,
        paramset = NA
    )

cat("\n>>> Writing output files...\n")
cat("Writing results to '", par$output_results, "'...\n", sep = "")
jsonlite::write_json(
    results,
    par$output_results,
    pretty = TRUE,
    null = "null",
    na = "null"
)

cat("\n>>> Done!\n")
