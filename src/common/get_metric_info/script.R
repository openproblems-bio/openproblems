library(tidyverse)
library(rlang)

## VIASH START
par <- list(
  input = "src",
  task_id = "label_projection",
  output = "output/metrics.json"
)
## VIASH END

ns_list <- processx::run(
  "viash",
  c("ns", "list", "-q", "metrics", "--src", paste("src", par$task_id, sep = "/")),
  wd = par$input
)
configs <- yaml::yaml.load(ns_list$stdout)

df <- map_df(configs, function(config) {
  if (length(config$functionality$status) > 0 && config$functionality$status == "disabled") return(NULL)
  info <- as_tibble(map_df(config$functionality$info$metrics, as.data.frame))
  info$config_path <- gsub(".*\\./", "", config$info$config)
  info$task_id <- par$task_id
  info$component_id <- config$functionality$name
  info$namespace <- config$functionality$namespace
  info$component_description <- config$functionality$description
  info$v1_url <- config$functionality$info$v1_url
  info$v1_commit <- config$functionality$info$v1_commit
  info
}) %>%
  select(metric_id, everything())

jsonlite::write_json(
  purrr::transpose(df),
  par$output,
  auto_unbox = TRUE,
  pretty = TRUE
)