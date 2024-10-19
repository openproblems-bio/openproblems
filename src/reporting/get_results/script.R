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
dir <- "work/c1/6660ea0cc6155d7e13fa341d16057b/_viash_par"
par <- list(
  input_scores = paste0(dir, "/input_scores_1/score_uns.yaml"),
  input_execution = paste0(dir, "/input_execution_1/trace.txt"),
  input_dataset_info = paste0(dir, "/input_dataset_info_1/output.json"),
  input_method_info = paste0(dir, "/input_method_info_1/output.json"),
  input_metric_info = paste0(dir, "/input_metric_info_1/output.json"),
  output_results = "output/results.json",
  output_metric_execution_info = "output/metric_execution_info.json"
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
      message("Encountered error while reading scores: ", e$message)
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

# read nxf log and process the task id
norm_methods <- "/log_cp10k|/log_cpm|/sqrt_cp10k|/sqrt_cpm|/l1_sqrt|/log_scran_pooling"
id_regex <- paste0("^.*:(.*)_process \\(([^\\.]*)(", norm_methods, ")?(.[^\\.]*)?\\.(.*)\\)$")

trace <- readr::read_tsv(par$input_execution) %>%
  mutate(
    id = name,
    process_id = stringr::str_extract(id, id_regex, 1L),
    dataset_id = stringr::str_extract(id, id_regex, 2L),
    normalization_id = gsub("^/", "", stringr::str_extract(id, id_regex, 3L)),
    grp4 = gsub("^\\.", "", stringr::str_extract(id, id_regex, 4L)),
    grp5 = stringr::str_extract(id, id_regex, 5L),
    submit = strptime(submit, "%Y-%m-%d %H:%M:%S"),
  ) %>%
  # detect whether entry is a metric or a method
  mutate(
    method_id = ifelse(is.na(grp4), grp5, grp4),
    metric_id = ifelse(is.na(grp4), grp4, grp5)
  ) %>%
  select(-grp4, -grp5) %>%
  filter(!is.na(method_id)) %>%
  # take last entry for each run
  arrange(desc(submit)) %>%
  group_by(name) %>%
  slice(1) %>%
  ungroup()

# parse values
execution_info <- trace %>%
  filter(process_id == method_id) %>% # only keep method entries
  rowwise() %>%
  transmute(
    dataset_id,
    normalization_id,
    method_id,
    resources = list(list(
      exit_code = parse_exit(exit),
      duration_sec = parse_duration(realtime),
      cpu_pct = parse_cpu(`%cpu`),
      peak_memory_mb = parse_size(peak_vmem),
      disk_read_mb = parse_size(rchar),
      disk_write_mb = parse_size(wchar)
    ))
  ) %>%
  ungroup()

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
metric_execution_info <- trace %>%
  filter(process_id == metric_id) %>% # only keep metric entries
  rowwise() %>%
  transmute(
    dataset_id,
    normalization_id,
    method_id,
    metric_id,
    resources = list(list(
      exit_code = parse_exit(exit),
      duration_sec = parse_duration(realtime),
      cpu_pct = parse_cpu(`%cpu`),
      peak_memory_mb = parse_size(peak_vmem),
      disk_read_mb = parse_size(rchar),
      disk_write_mb = parse_size(wchar)
    ))
  ) %>%
  ungroup()

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
