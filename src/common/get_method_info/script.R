library(tidyverse)
library(rlang)

## VIASH START
par <- list(
  input = "src/label_projection",
  output = "resources/label_projection/output/method_info.yaml"
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
  info$id <- config$functionality$name
  info$namespace <- config$functionality$namespace
  info$description <- config$functionality$description
  info
}) %>%
  select(id, type, label, everything())

yaml::write_yaml(purrr::transpose(df), par$output)

