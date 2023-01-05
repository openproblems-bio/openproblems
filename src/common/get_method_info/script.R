library(tidyverse)
library(rlang)

## VIASH START
par <- list(
  input = "src",
  task_id = "label_projection",
  output = "output/method_info.json"
)
## VIASH END

ns_list <- processx::run(
  "viash",
  c("ns", "list", "-q", "methods", "--src", paste("src", par$task_id, sep = "/")),
  wd = par$input
)
configs <- yaml::yaml.load(ns_list$stdout)

df <- map_df(configs, function(config) {
  if (length(config$functionality$status) > 0 && config$functionality$status == "disabled") return(NULL)
  info <- as_tibble(config$functionality$info)
  info$config_path <- gsub(".*\\./", "", config$info$config)
  info$id <- config$functionality$name
  info$namespace <- config$functionality$namespace
  info$description <- config$functionality$description
  info
}) %>%
  select(id, type, label, everything())

jsonlite::write_json(
  purrr::transpose(df),
  par$output,
  auto_unbox = TRUE,
  pretty = TRUE
)