library(tidyverse)
library(rlang)

## VIASH START
par <- list(
  input = "src/label_projection",
  output = "resources/label_projection/output/metric_info.yaml"
)
## VIASH END

ns_list <- processx::run(
  "viash",
  c("ns", "list", "-q", "metrics", "--src", "."),
  wd = par$input
)
configs <- yaml::yaml.load(ns_list$stdout)

df <- map_df(configs, function(config) {
  info <- as_tibble(map_df(config$functionality$info$metrics, as.data.frame))
  info$config_path <- gsub(".*\\./", "", config$info$config)
  info$id <- config$functionality$name
  info$namespace <- config$functionality$namespace
  info$description <- config$functionality$description
  info
})

yaml::write_yaml(purrr::transpose(df), par$output)

