requireNamespace("jsonlite", quietly = TRUE)
requireNamespace("yaml", quietly = TRUE)
requireNamespace("dynutils", quietly = TRUE)
requireNamespace("readr", quietly = TRUE)
requireNamespace("lubridate", quietly = TRUE)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)
library(rlang, warn.conflicts = FALSE)

## VIASH START
# raw_dir <- "resources_test/openproblems/task_results_v3/raw"
# processed_dir <- "resources_test/openproblems/task_results_v3/processed"
# raw_dir <- "/home/rcannood/workspace/openproblems-bio/task_perturbation_prediction/resources/results/run_2024-10-31_06-14-14"
# processed_dir <- "/home/rcannood/workspace/openproblems-bio/website/results/perturbation_prediction/data"
raw_dir <- "/home/rcannood/workspace/openproblems-bio/task_batch_integration/resources/results/run_2024-11-20_12-47-03"
processed_dir <- "/home/rcannood/workspace/openproblems-bio/website/results/batch_integration/data"

par <- list(
  # inputs
  input_scores = paste0(raw_dir, "/score_uns.yaml"),
  input_execution = paste0(raw_dir, "/trace.txt"),
  input_dataset_info = paste0(processed_dir, "/dataset_info.json"),
  input_method_info = paste0(processed_dir, "/method_info.json"),
  input_method_configs = paste0(raw_dir, "/method_configs.yaml"),
  input_metric_info = paste0(processed_dir, "/metric_info.json"),
  # outputs
  output_results = paste0(processed_dir, "/results.json"),
  output_metric_execution_info = paste0(processed_dir, "/metric_execution_info.json")
)
## VIASH END

# --- helper functions ---------------------------------------------------------
cat("Loading helper functions\n")
parse_exit <- function(x) {
  if (is.na(x) || x == "-") {
    NA_integer_
  } else {
    as.integer(x)
  }
}
parse_duration <- function(x) {
  if (is.na(x) || x == "-") {
    NA_real_
  } else {
    as.numeric(lubridate::duration(toupper(x)))
  }
}
parse_cpu <- function(x) {
  if (is.na(x) || x == "-") {
    NA_real_
  } else {
    as.numeric(gsub(" *%", "", x))
  }
}
parse_size <- function(x) {
  out <-
    if (is.na(x) || x == "-") {
      NA_integer_
    } else if (grepl("TB", x)) {
      as.numeric(gsub(" *TB", "", x)) * 1024 * 1024
    } else if (grepl("GB", x)) {
      as.numeric(gsub(" *GB", "", x)) * 1024
    } else if (grepl("MB", x)) {
      as.numeric(gsub(" *MB", "", x))
    } else if (grepl("KB", x)) {
      as.numeric(gsub(" *KB", "", x)) / 1024
    } else if (grepl("B", x)) {
      as.numeric(gsub(" *B", "", x)) / 1024 / 1024
    } else {
      NA_integer_
    }
  as.integer(ceiling(out))
}

# --- read input files ---------------------------------------------------------
cat("Reading input files\n")
# read scores
raw_scores <-
  yaml::yaml.load_file(par$input_scores) %>%
  map_df(function(x) {
    tryCatch({
      as_tibble(as.data.frame(
        x[c("dataset_id", "method_id", "metric_ids", "metric_values")]
      ))
    }, error = function(e) {
      message("Encountered error while reading scores.\n  Error: ", e$message, "\n  Data: ", paste(paste0(names(x), "=", x), collapse = ", "))
      NULL
    })
  })

# read metric info
dataset_info <- jsonlite::read_json(par$input_dataset_info, simplifyVector = TRUE)
method_info <- jsonlite::read_json(par$input_method_info, simplifyVector = TRUE)
metric_info <- jsonlite::read_json(par$input_metric_info, simplifyVector = TRUE)

# --- process scores and execution info ----------------------------------------
cat("Processing scores and execution info\n")
scale_scores <- function(values, is_control, maximize) {
  control_values <- values[is_control & !is.na(values)]
  if (length(control_values) < 2) {
    return(NA_real_)
  }

  min_control_value <- min(control_values)
  max_control_value <- max(control_values)

  if (min_control_value == max_control_value) {
    return(NA_real_)
  }

  scaled <- (values - min_control_value) / (max_control_value - min_control_value)

  if (maximize) {
    scaled
  } else {
    1 - scaled
  }
}
aggregate_scores <- function(scaled_score) {
  mean(pmin(1, pmax(0, scaled_score)) %|% 0)
}
scores <- raw_scores %>%
  complete(
    dataset_id,
    method_id,
    metric_ids,
    fill = list(metric_values = NA_real_)
  ) %>%
  left_join(method_info %>% select(method_id, is_baseline), by = "method_id") %>%
  left_join(metric_info %>% select(metric_ids = metric_id, maximize), by = "metric_ids") %>%
  group_by(metric_ids, dataset_id) %>%
  mutate(scaled_score = scale_scores(metric_values, is_baseline, maximize[[1]]) %|% 0) %>%
  group_by(dataset_id, method_id) %>%
  summarise(
    metric_values = list(as.list(setNames(metric_values, metric_ids))),
    scaled_scores = list(as.list(setNames(scaled_score, metric_ids))),
    mean_score = aggregate_scores(scaled_score),
    .groups = "drop"
  )


# read execution info
# -> only keep the last execution of each process
input_execution <- readr::read_tsv(par$input_execution) |>
  group_by(name) |>
  mutate(num_runs = n()) |>
  slice(which.max(submit)) |>
  ungroup()

method_lookup <- map_dfr(method_info$method_id, function(method_id) {
  regex <- paste0("(.*:", method_id, ":[^ ]*)")
  name <-
    input_execution$name[grepl(regex, input_execution$name)] |>
    unique()
  name_ <- name[!grepl(":publishStatesProc", name)]
  tibble(method_id = method_id, name = name_)
})
dataset_lookup <- map_dfr(dataset_info$dataset_id, function(dataset_id) {
  regex <- paste0(".*[(.](", dataset_id, ")[)./].*")
  name <-
    input_execution$name[grepl(regex, input_execution$name)] |>
    unique()
  tibble(dataset_id = dataset_id, name = name)
})

# parse values
execution_info_ind <- input_execution |>
  left_join(method_lookup, by = "name") |>
  left_join(dataset_lookup, by = "name") |>
  filter(!is.na(method_id)) %>%
  rowwise() |>
  mutate(
    process_id = gsub(" .*", "", name),
    submit = strptime(submit, "%Y-%m-%d %H:%M:%S"),
    exit_code = parse_exit(exit),
    duration_sec = parse_duration(realtime),
    cpu_pct = parse_cpu(`%cpu`),
    peak_memory_mb = parse_size(peak_vmem),
    disk_read_mb = parse_size(rchar),
    disk_write_mb = parse_size(wchar)
  ) |>
  ungroup()

execution_info <- execution_info_ind |>
  group_by(dataset_id, method_id) |>
  summarise(
    resources = list(list(
      submit = min(submit),
      exit_code = max(exit_code),
      duration_sec = sum(duration_sec),
      cpu_pct = sum(cpu_pct * duration_sec) / sum(duration_sec),
      peak_memory_mb = max(peak_memory_mb),
      disk_read_mb = sum(disk_read_mb),
      disk_write_mb = sum(disk_write_mb)
    )),
    .groups = "drop"
  )

# combine scores with execution info
# fill up missing entries with NAs and 0s
metric_ids <- unique(raw_scores$metric_ids)
rep_names <- function(val) {
  setNames(
    as.list(rep(val, length(metric_ids))),
    metric_ids
  )
}
out <- full_join(
  scores,
  execution_info,
  by = c("method_id", "dataset_id")
) %>%
  rowwise() %>%
  mutate(
    task_id = par$task_id,
    metric_values = list(metric_values %||% rep_names(NA_real_)),
    scaled_scores = list(scaled_scores %||% rep_names(0)),
    mean_score = mean_score %|% 0,
  ) %>%
  ungroup()


# --- process metric execution info --------------------------------------------
cat("Processing metric execution info\n")

# manually add component id to metric info
metric_info$component_name <- metric_info$component_name %||% rep(NA_character_, nrow(metric_info)) %|%
  gsub(".*/([^/]*)/config\\.vsh\\.yaml", "\\1", metric_info$implementation_url)

metric_lookup2 <- pmap_dfr(metric_info, function(metric_id, component_name, ...) {
  regex <- paste0("(.*:", component_name, ":[^ ]*)")
  name <-
    input_execution$name[grepl(regex, input_execution$name)] |>
    unique()
  name_ <- name[!grepl(":publishStatesProc", name)]
  tibble(metric_id = metric_id, component_name = component_name, name = name_)
})
dataset_lookup2 <- map_dfr(dataset_info$dataset_id, function(dataset_id) {
  regex <- paste0(".*[(.](", dataset_id, ")[)./].*")
  name <-
    input_execution$name[grepl(regex, input_execution$name)] |>
    unique()
  tibble(dataset_id = dataset_id, name = name)
})
method_lookup2 <- map_dfr(method_info$method_id, function(method_id) {
  regex <- paste0(".*[(.](", method_id, ")[)./].*")
  name <-
    input_execution$name[grepl(regex, input_execution$name)] |>
    unique()
  tibble(method_id = method_id, name = name)
})

metric_execution_info_ind <- input_execution |>
  left_join(metric_lookup2, by = "name") |>
  left_join(dataset_lookup2, by = "name") |>
  left_join(method_lookup2, by = "name") |>
  filter(!is.na(metric_id)) %>%
  rowwise() |>
  mutate(
    process_id = gsub(" .*", "", name),
    submit = strptime(submit, "%Y-%m-%d %H:%M:%S"),
    exit_code = parse_exit(exit),
    duration_sec = parse_duration(realtime),
    cpu_pct = parse_cpu(`%cpu`),
    peak_memory_mb = parse_size(peak_vmem),
    disk_read_mb = parse_size(rchar),
    disk_write_mb = parse_size(wchar)
  ) |>
  ungroup()

metric_execution_info <- metric_execution_info_ind |>
  group_by(dataset_id, method_id, metric_component_name = component_name) |>
  summarise(
    resources = list(list(
      submit = min(submit),
      exit_code = max(exit_code),
      duration_sec = sum(duration_sec),
      cpu_pct = sum(cpu_pct * duration_sec) / sum(duration_sec),
      peak_memory_mb = max(peak_memory_mb),
      disk_read_mb = sum(disk_read_mb),
      disk_write_mb = sum(disk_write_mb)
    )),
    .groups = "drop"
  )


# --- write output files -------------------------------------------------------
cat("Writing output files\n")
# write output files
jsonlite::write_json(
  purrr::transpose(out),
  par$output_results,
  auto_unbox = TRUE,
  pretty = TRUE
)
jsonlite::write_json(
  purrr::transpose(metric_execution_info),
  par$output_metric_execution_info,
  auto_unbox = TRUE,
  pretty = TRUE
)
