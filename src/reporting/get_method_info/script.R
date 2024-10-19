requireNamespace("jsonlite", quietly = TRUE)
requireNamespace("yaml", quietly = TRUE)
library(purrr, warn.conflicts = FALSE)
library(rlang, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input = "resources_test/openproblems/task_results_v3/raw/method_configs.yaml",
  output = "resources_test/openproblems/task_results_v3/processed/method_info.json"
)
## VIASH END

configs <- yaml::yaml.load_file(par$input)

outputs <- map(configs, function(config) {
  if (length(config$functionality$status) > 0 && config$functionality$status == "disabled") {
    return(NULL)
  }

  # prep for viash 0.9.0
  build_info <- config$build_info %||% config$info
  if ("functionality" %in% names(config)) {
    config[names(config$functionality)] <- config$functionality
    config[["functionality"]] <- NULL
  }

  info <- config$info

  # add extra info
  info$config_path <- gsub(".*/src/", "src/", build_info$config)
  info$task_id <- gsub("/.*", "", config$namespace)
  info$id <- config$name
  info$namespace <- config$namespace
  info$label <- config$label %||% info$label
  info$summary <- config$summary %||% info$summary
  info$description <- config$description %||% info$description
  info$commit_sha <- build_info$git_commit %||% "missing-sha"
  info$code_version <- "missing-version"
  info$implementation_url <- paste0(
    build_info$git_remote, "/blob/",
    build_info$git_commit, "/",
    info$config_path
  )
  info$type_info <- NULL

  # ↑ this could be used as the new format

  # construct v1 format
  out <- list(
    task_id = info$task_id,
    method_id = info$id,
    method_name = info$label,
    method_summary = info$summary,
    method_description = info$description,
    is_baseline = grepl("control", info$type),
    paper_reference = info$reference %||% NA_character_,
    code_url = info$repository_url %||% NA_character_,
    implementation_url = info$implementation_url %||% NA_character_,
    code_version = NA_character_,
    commit_sha = info$commit_sha
  )

  # show warning when certain data is missing and return null?
  for (n in names(out)) {
    if (is.null(out[[n]])) {
      out_as_str <- jsonlite::toJSON(out, auto_unbox = TRUE, pretty = TRUE)
      stop("missing value for value '", n, "' in ", out_as_str)
    }
  }

  # return output
  out
})

jsonlite::write_json(
  outputs,
  par$output,
  auto_unbox = TRUE,
  pretty = TRUE
)