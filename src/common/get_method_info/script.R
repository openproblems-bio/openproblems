library(tidyverse)
library(rlang)

## VIASH START
par <- list(
  input = ".",
  task_id = "label_projection",
  output = "output/method_info.json"
)
## VIASH END

ns_list <- processx::run(
  "viash",
  c("ns", "list", "-q", "methods", "--src", paste("src/tasks", par$task_id, sep = "/")),
  wd = par$input
)
configs <- yaml::yaml.load(ns_list$stdout)

out <- map(configs, function(config) {
  if (length(config$functionality$status) > 0 && config$functionality$status == "disabled") return(NULL)
  info <- config$functionality$info

  # add extra info
  info$config_path <- gsub(".*\\./", "", config$info$config)
  info$task_id <- par$task_id
  info$method_id <- config$functionality$name
  info$namespace <- config$functionality$namespace
  info$is_baseline <- grepl("control", info$type)

  # rename fields to v1 format
  info$method_name <- info$pretty_name
  info$pretty_name <- NULL
  info$method_summary <- info$description
  info$description <- NULL
  info$paper_reference <- info$reference
  info$reference <- NULL
  info$code_url <- info$repository_url
  info$repository_url <- NULL

  # todo: show warning when certain data is missing and return null?

  # return output
  info
})

jsonlite::write_json(
  out,
  par$output,
  auto_unbox = TRUE,
  pretty = TRUE
)