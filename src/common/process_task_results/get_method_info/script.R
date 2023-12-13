library(purrr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(rlang, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input = "output/temp/method_configs.yaml",
  task_id = "label_projection",
  output = "output/test/method_info.json"
)
## VIASH END

configs <- yaml::yaml.load_file(par$input)

out <- map(configs, function(config) {
  if (length(config$functionality$status) > 0 && config$functionality$status == "disabled") return(NULL)
  info <- config$functionality$info

  # add extra info
  info$config_path <- gsub(".*\\./", "", config$info$config)
  info$task_id <- par$task_id
  info$method_id <- config$functionality$name
  info$namespace <- config$functionality$namespace
  info$is_baseline <- grepl("control", info$type)
  info$commit_sha <- config$info$git_commit %||% "missing-sha"
  info$code_version <- "missing-version"

  # rename fields to v1 format
  info$method_name <- info$label
  info$label <- NULL
  info$method_summary <- info$summary
  info$summary <- NULL
  info$method_description <- info$description
  info$description <- NULL
  info$paper_reference <- info$reference
  info$reference <- NULL
  info$code_url <- info$repository_url
  info$repository_url <- NULL
  info$v1.path <- info$v1$path
  info$v1$path <- NULL
  info$v1.commit <- info$v1$commit
  info$v1$commit <- NULL
  info$v1 <- NULL
  info$type_info.label <- info$type_info$label
  info$type_info$label <- NULL
  info$type_info.summary <- info$type_info$summary
  info$type_info$summary <- NULL
  info$type_info.description <- info$type_info$description
  info$type_info$description <- NULL
  info$type_info <- NULL
  if (length(info$variants) > 0) {
    info$variants <- NULL
  }



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