requireNamespace("jsonlite", quietly = TRUE)
requireNamespace("yaml", quietly = TRUE)
library(purrr, warn.conflicts = FALSE)
library(rlang, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input = "resources_test/openproblems/task_results_v3/raw/metric_configs.yaml",
  output = "output/metric_info.json"
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

  map(
    config$info$metrics,
    function(info) {
      # add extra info
      info$config_path <- gsub(".*/src/", "src/", build_info$config)
      info$task_id <- gsub("/.*", "", config$namespace)
      info$id <- info$name
      info$name <- NULL
      info$component_id <- config$name
      info$namespace <- config$namespace
      info$commit_sha <- build_info$git_commit %||% "missing-sha"
      info$code_version <- "missing-version"
      info$implementation_url <- paste0(
        build_info$git_remote, "/blob/",
        build_info$git_commit, "/",
        info$config_path
      )

      # ↑ this could be used as the new format

      # construct v1 format
      out <- list(
        task_id = info$task_id,
        metric_id = info$id,
        metric_name = info$label,
        metric_summary = info$summary,
        metric_description = info$description,
        paper_reference = info$reference %||% NA_character_,
        implementation_url = info$implementation_url %||% NA_character_,
        code_version = NA_character_,
        commit_sha = info$commit_sha,
        maximize = info$maximize
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
    }
  )
})

outputs <- unlist(outputs, recursive = FALSE)

jsonlite::write_json(
  outputs,
  par$output,
  auto_unbox = TRUE,
  pretty = TRUE
)