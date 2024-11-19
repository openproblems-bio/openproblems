requireNamespace("jsonlite", quietly = TRUE)
requireNamespace("yaml", quietly = TRUE)
library(purrr, warn.conflicts = FALSE)
library(rlang, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input = "test.yaml",
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
  info$comp_path <- gsub(".*/src/", "src/", build_info$config) %>% gsub("/config.vsh.yaml", "", .)
  info$task_id <- gsub("/.*", "", config$namespace)
  info$id <- config$name
  info$namespace <- config$namespace
  info$label <- config$label %||% info$label
  info$summary <- config$summary %||% info$summary
  info$description <- config$description %||% info$description
  info$commit_sha <- build_info$git_commit %||% "missing-sha"
  info$code_version <- config$version
  info$code_url <- config$links$repository
  info$documentation_url <- config$links$documentation
  # Check if the method has a docker container to create an image url. If it does not have a docker it will be a nextflow component consisting of different components that will have a docker image.
  engines <- config$engines
  has_docker <- any(map_lgl(engines, ~ .x$type == "docker"))
  if (has_docker) {
    info$image <- paste0(
      "https://",
      config$links$docker_registry, "/",
      config$package_config$organization, "/",
      config$package_config$name, "/",
      gsub("src/", "", info$comp_path),
      ":",
      info$code_version
    )
  }  else {
    info$image <- paste0(
      "https://github.com/orgs/openproblems-bio/packages?repo_name=",
      config$package_config$name,
      "&q=",
      gsub("src/", "", info$comp_path)
    )
  }
  info$implementation_url <- paste0(
    build_info$git_remote, "/blob/",
    build_info$git_commit, "/",
    info$comp_path
  )
  info$type_info <- NULL

  # Flatten references
  if (!is.null(info$reference) && info$reference != "") {
    imap(info$reference, function(value, key) {
      info[[paste0("reference_", key)]] <- value
    })
  }
  info$reference <- NULL


  # ↑ this could be used as the new format

  # construct v1 format
  out <- list(
    task_id = info$task_id,
    method_id = info$id,
    method_name = info$label,
    method_summary = info$summary,
    method_description = info$description,
    is_baseline = grepl("control", info$type),
    reference_doi = info$reference_doi %||% NA_character_,
    reference_bibtex = info$reference_bibtex %||% NA_character_,
    code_url = info$code_url %||% NA_character_,
    documentation_url = info$documentation_url %||% NA_character_,
    image = info$image %||% NA_character_,
    implementation_url = info$implementation_url %||% NA_character_,
    code_version = info$code_version %||% NA_character_,
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