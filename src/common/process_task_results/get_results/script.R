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
par <- list(
  input_scores = "output/temp/score_uns.yaml",
  input_execution = "output/temp/trace.txt",
  output = "output/results.json"
)
## VIASH END

# read scores
raw_scores <- yaml::yaml.load_file(par$input_scores) %>%
  map_df(function(x) {
    as_tibble(as.data.frame(
      x[c("dataset_id", "method_id", "metric_ids", "metric_values")]
    ))
  })

# scale scores
scores <- raw_scores %>%
  complete(
    dataset_id,
    method_id,
    metric_ids,
    fill = list(metric_values = NA_real_)
  ) %>%
  group_by(metric_ids, dataset_id) %>%
  mutate(
    scaled_score = dynutils::scale_minmax(metric_values) %|% 0
  ) %>%
  group_by(dataset_id, method_id) %>%
  summarise(
    metric_values = list(as.list(setNames(metric_values, metric_ids))),
    scaled_scores = list(as.list(setNames(scaled_score, metric_ids))),
    mean_score = mean(scaled_score),
    .groups = "drop"
  )

# read nxf log and process the task id
id_regex <- "^.*:(.*)_process \\((.*)/([^\\.]*)\\.(.*)\\)$"

trace <- readr::read_tsv(par$input_execution) %>%
  mutate(
    id = name,
    process_id = stringr::str_extract(id, id_regex, 1L),
    dataset_id = stringr::str_extract(id, id_regex, 2L),
    normalization_id = stringr::str_extract(id, id_regex, 3L),
    method_id = stringr::str_extract(id, id_regex, 4L),
  ) %>%
  filter(process_id == method_id)

# parse strings into numbers
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

execution_info <- trace %>%
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


jsonlite::write_json(
  purrr::transpose(out),
  par$output,
  auto_unbox = TRUE,
  pretty = TRUE
)
