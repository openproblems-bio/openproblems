library(tidyverse)
library(rlang)

## VIASH START
par <- list(
  input = "src/denoising",
  output = "temp/denoising_metrics.yaml"
)
## VIASH END

ns_list <- processx::run(
  "viash",
  c("ns", "list", "-q", "metrics", "--src", "."),
  wd = par$input
)
configs <- yaml::yaml.load(ns_list$stdout)

df <- map_df(configs, function(config) {
  if (length(config$functionality$status) > 0 && config$functionality$status == "disabled") return(NULL)
  info <- as_tibble(map_df(config$functionality$info$metrics, as.data.frame))
  info$config_path <- gsub(".*\\./", "", config$info$config)
  info$component_id <- config$functionality$name
  info$namespace <- config$functionality$namespace
  info$component_description <- config$functionality$description
  info$v1_url <- config$functionality$info$v1_url
  info$v1_commit <- config$functionality$info$v1_commit
  info
}) %>%
  select(id, everything())

yaml::write_yaml(purrr::transpose(df), par$output)