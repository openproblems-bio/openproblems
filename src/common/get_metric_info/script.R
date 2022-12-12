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
  colnames(info) <- paste0("metrics_", colnames(info))
  info$config_path <- gsub(".*\\./", "", config$info$config)
  info$id <- config$functionality$name # overwrites metrics id first add the Id to metrics _id ?
  info$namespace <- config$functionality$namespace
  info$description <- config$functionality$description # same as id 
  info$v1_url <- config$functionality$info$v1_url
  info$v1_commit <- config$functionality$info$v1_commit
  info
}) %>%
  select(id, everything())

yaml::write_yaml(purrr::transpose(df), par$output)