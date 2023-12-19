library(purrr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(rlang, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input = "resources_test/common/task_metadata/dataset_info.yaml",
  output = "output/metric_info.json"
)
## VIASH END

datasets <- yaml::yaml.load_file(par$input)

df <- map_df(datasets, function(dataset) {
  info <- as_tibble(map(dataset, as.data.frame))
}) %>%
  rename(
    data_url = dataset_url,
    data_reference = dataset_reference
  )


jsonlite::write_json(
  purrr::transpose(df),
  par$output,
  auto_unbox = TRUE,
  pretty = TRUE
)