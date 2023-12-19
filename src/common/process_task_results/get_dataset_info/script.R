requireNamespace("jsonlite", quietly = TRUE)
requireNamespace("yaml", quietly = TRUE)

## VIASH START
par <- list(
  input = "resources_test/common/task_metadata/dataset_info.yaml",
  output = "output/dataset_info.json"
)
## VIASH END

datasets <- yaml::yaml.load_file(par$input)

# transform into format expected by website
datasets_formatted <- lapply(datasets, function(dataset) {
  dataset$data_url <- dataset$dataset_url
  dataset$data_reference <- dataset$dataset_reference
  dataset
})

jsonlite::write_json(
  datasets_formatted,
  par$output,
  auto_unbox = TRUE,
  pretty = TRUE
)