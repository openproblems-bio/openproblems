library(purrr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(rlang, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input = "output/temp/metric_configs.yaml",
  task_id = "batch_integration",
  output = "output/metric_info.json"
)
## VIASH END

configs <- yaml::yaml.load_file(par$input)

df <- map_df(configs, function(config) {
  if (length(config$functionality$status) > 0 && config$functionality$status == "disabled") return(NULL)
  info <- as_tibble(map_df(config$functionality$info$metrics, as.data.frame))
  info$config_path <- gsub(".*\\./", "", config$info$config)
  info$task_id <- par$task_id
  info$component_id <- config$functionality$name
  info$namespace <- config$functionality$namespace
  info
}) %>%
  rename(
    metric_id = name,
    metric_name = label,
    metric_summary = description,
    paper_reference = reference,
  )

jsonlite::write_json(
  purrr::transpose(df),
  par$output,
  auto_unbox = TRUE,
  pretty = TRUE
)