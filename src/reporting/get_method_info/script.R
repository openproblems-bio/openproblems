## VIASH START
par <- list(
  input = "resources_test/openproblems/task_results_v4/processed/task_info.json",
  output = "resources_test/openproblems/task_results_v4/processed/method_info.json"
)
## VIASH END

################################################################################
#                               FUNCTIONS
################################################################################

get_implementation_url <- function(config) {
  paste0(
    config$build_info$git_remote,
    "/blob/",
    config$build_info$git_commit,
    "/",
    config$build_info$config |>
      stringr::str_replace(".*/src/", "src/") |>
      stringr::str_remove("/config.vsh.yaml")
  )
}

get_container_image <- function(config) {
  # Check if the method has a docker container to create an image url.
  # If it does not have a docker it will be a nextflow component consisting of
  # different components that will have a docker image.
  engines <- config$engines
  has_docker <- any(purrr::map_lgl(engines, ~ .x$type == "docker"))
  if (has_docker) {
    paste0(
      "https://",
      config$links$docker_registry, "/",
      config$package_config$organization, "/",
      config$package_config$name, "/",
      config$build_info$config |>
        stringr::str_remove(".*/src/") |>
        stringr::str_remove("/config.vsh.yaml"),
      ":",
      config$version
    )
  }  else {
    paste0(
      "https://github.com/orgs/openproblems-bio/packages?repo_name=",
      config$package_config$name,
      "&q=",
      config$build_info$config |>
        stringr::str_remove(".*/src/") |>
        stringr::str_remove("/config.vsh.yaml")
    )
  }
}

get_references <- function(config) {
  if (!is.null(config$references)) {
    list(
      doi = config$references$doi %||% character(0),
      bibtex = config$references$bibtex %||% character(0)
    )
  } else {
    list(doi = character(0), bibtex = character(0))
  }
}

get_additional_info <- function(config) {
  # Fields that are stored elsewhere and we don't want to save here
  exclude <- c("type", "type_info")

  config$info[setdiff(names(config$info), exclude)] |>
    purrr::map(recurse_unbox)
}

recurse_unbox <- function(x) {
  if (is.list(x)) {
    purrr::map(x, recurse_unbox)
  } else if (length(x) == 1) {
    jsonlite::unbox(x)
  } else {
    x
  }
}

################################################################################
#                              MAIN SCRIPT
################################################################################

cat("====== Get method info ======\n")

`%||%` <- rlang::`%||%`

cat("\n>>> Reading input files...\n")
cat("Reading method info from '", par$input, "'...\n", sep = "")
method_configs <- yaml::yaml.load_file(par$input)

cat("\n>>> Processing ", length(method_configs), " method configs...\n", sep = "")
method_info_json <- purrr::map(method_configs, function(.config) {
  if (.config$status == "disabled") {
    cat("Skipping disabled method '", .config$name, "'\n", sep = "")
    return(NULL)
  } else {
    cat("Processing method '", .config$name, "'\n", sep = "")
  }

  list(
    name = jsonlite::unbox(.config$name),
    label = jsonlite::unbox(.config$label),
    commit = jsonlite::unbox(.config$build_info$git_commit %||% "missing-sha"),
    summary = .config$summary |>
      stringr::str_trim() |>
      stringr::str_remove_all('(^"|"$|^\'|\'$)') |>
      jsonlite::unbox(),
    description = .config$description |>
      stringr::str_trim() |>
      stringr::str_remove_all('(^"|"$|^\'|\'$)') |>
      jsonlite::unbox(),
    type = jsonlite::unbox(.config$info$type),
    link_code = jsonlite::unbox(.config$links$repository),
    link_documentation = jsonlite::unbox(.config$links$documentation),
    link_implementation = jsonlite::unbox(get_implementation_url(.config)),
    link_container_image = jsonlite::unbox(get_container_image(.config)),
    references = get_references(.config),
    additional_info = get_additional_info(.config),
    version = jsonlite::unbox(.config$version)
  )
})

cat("\n>>> Writing output files...\n")
cat("Writing task info to '", par$output, "'...\n", sep = "")
jsonlite::write_json(
  method_info_json,
  par$output,
  pretty = TRUE,
  null = "null"
)

cat("\n>>> Validating output against schema...\n")
ajv_args <- paste(
  "validate",
  "--spec draft2020",
  "-s", file.path(meta$resources_dir, "schemas", "method_info_schema.json"),
  "-r", file.path(meta$resources_dir, "schemas", "references_schema.json"),
  "-d", par$output
)

cat("Running validation command:", "ajv", ajv_args, "\n")
cat("Output:\n")
validation_result <- system2("ajv", ajv_args)

if (validation_result == 0) {
  cat("JSON validation passed successfully!\n")
} else {
  cat("JSON validation failed!\n")
  stop("Output JSON does not conform to schema")
}

cat("\n>>> Done!\n")
