requireNamespace("jsonlite", quietly = TRUE)
requireNamespace("yaml", quietly = TRUE)
library(purrr, warn.conflicts = FALSE)
library(rlang, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input = "output/temp/metric_configs.yaml",
  output = "output/metric_info.json"
)
## VIASH END

configs <- yaml::yaml.load_file(par$input)

outputs <- map(configs, function(config) {
  if (length(config$functionality$status) > 0 && config$functionality$status == "disabled") {
    return(NULL)
  }

  map(
    config$functionality$info$metrics,
    function(info) {
      # add extra info
      info$config_path <- gsub(".*openproblems-v2/src/", "src/", config$info$config)
      info$task_id <- gsub("/.*", "", config$functionality$namespace)
      info$id <- info$name
      info$component_id <- config$functionality$name
      info$namespace <- config$functionality$namespace
      info$commit_sha <- config$info$git_commit %||% "missing-sha"
      info$code_version <- "missing-version"
      info$implementation_url <- paste0(
        "https://github.com/openproblems-bio/openproblems-v2/tree/",
        info$commit_sha, "/",
        info$config_path
      )

      # ↑ this could be used as the new format

      # construct v1 format
      out <- list(
        task_id = info$task_id,
        metric_id = info$id,
        metric_name = info$label,
        metric_summary = info$description,
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