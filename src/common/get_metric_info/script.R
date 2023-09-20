library(purrr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(rlang, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input = ".",
  task_id = "batch_integration",
  output = "output/metric_info.json"
)
## VIASH END

ns_list <- processx::run(
  "viash",
  c("ns", "list", "-q", "metrics", "--src", paste("src/tasks", par$task_id, sep = "/")),
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