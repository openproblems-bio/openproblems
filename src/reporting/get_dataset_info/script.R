requireNamespace("jsonlite", quietly = TRUE)
requireNamespace("yaml", quietly = TRUE)
library(purrr, warn.conflicts = FALSE)
library(rlang, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input = "output/label_projection/dataset_uns.yaml",
  output = "output/dataset_info.json"
)
## VIASH END

datasets <- yaml::yaml.load_file(par$input)

# transform into format expected by website
outputs <- map(datasets, function(dataset) {
  # ↑ the 'dataset' object could be used as the new format

  # TODO: it'd be nice if the s3 path was also included in the dataset info

  # construct v1 format
  out <- list(
    "dataset_id" = dataset$dataset_id,
    "dataset_name" = dataset$dataset_name,
    "dataset_summary" = dataset$dataset_summary,
    "dataset_description" = dataset$dataset_description %||% NA_character_,
    "data_reference" = dataset$dataset_reference %||% NA_character_,
    "data_url" = dataset$dataset_url %||% NA_character_,
    "date_created" = dataset$date_created %||% NA_character_,
    "file_size" = dataset$file_size %||% NA_character_
  )

  if (!is.null(dataset[["common_dataset_id"]])) {
    out[["common_dataset_id"]] <- dataset[["common_dataset_id"]]
  }

  # show warning when certain data is missing and return null?
  for (n in names(out)) {
    if (is.null(out[[n]])) {
      out_as_str <- jsonlite::toJSON(out, auto_unbox = TRUE, pretty = TRUE)
      stop("missing value for value '", n, "' in ", out_as_str)
    }
  }

  out
})

jsonlite::write_json(
  outputs,
  par$output,
  auto_unbox = TRUE,
  pretty = TRUE
)
