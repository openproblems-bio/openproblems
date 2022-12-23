library(tidyverse)
library(rlang)

## VIASH START
par <- list(
  input = "src/label_projection",
  output = "resources_test/common/method_info.json"
)
## VIASH END

ns_list <- processx::run(
  "viash",
  c("ns", "list", "-q", "methods", "--src", "."),
  wd = par$input
)
configs <- yaml::yaml.load(ns_list$stdout)

df <- map_df(configs, function(config) {
  if (length(config$functionality$status) > 0 && config$functionality$status == "disabled") return(NULL)
  info <- as_tibble(config$functionality$info)
  info$config_path <- gsub(".*\\./", "", config$info$config)
  info$task_id <- par$query
  info$method_id <- config$functionality$name
  info$namespace <- config$functionality$namespace
  info$description <- config$functionality$description
  info$is_baseline <- FALSE
  if (grepl("control", info$type)) {
    info$is_baseline <- TRUE
  }
  info
}) %>%
  select(method_id, type, method_name, everything())

jsonlite::write_json(
  purrr::transpose(df),
  par$output,
  auto_unbox = TRUE,
  pretty = TRUE
)