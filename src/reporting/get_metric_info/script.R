## VIASH START
par <- list(
  input = "resources_test/openproblems/task_results_v4/raw/metric_configs.yaml",
  output = "resources_test/openproblems/task_results_v4/processed/metric_info.json"
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
  } else {
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

get_additional_info <- function(info, exclude) {
  info[setdiff(names(info), exclude)] |>
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

cat("====== Get metric info ======\n")

`%||%` <- rlang::`%||%`

cat("\n>>> Reading input files...\n")
cat("Reading metric info from '", par$input, "'...\n", sep = "")
metric_configs <- yaml::yaml.load_file(par$input)

cat(
  "\n>>> Processing ", length(metric_configs), " metric configs...\n",
  sep = ""
)
metric_info_json <- purrr::map(metric_configs, function(.config) {
  if (.config$status == "disabled") {
    cat("Skipping disabled metric component '", .config$name, "'\n", sep = "")
    return(NULL)
  } else {
    cat("Processing metric component '", .config$name, "'\n", sep = "")
  }

  purrr::map(.config$info$metrics, function(.metric) {
    list(
      name = jsonlite::unbox(.metric$name),
      label = jsonlite::unbox(.metric$label),
      commit = jsonlite::unbox(
        .config$build_info$git_commit %||% "missing-sha"
      ),
      summary = jsonlite::unbox(.metric$summary),
      description = jsonlite::unbox(.metric$description),
      maximize = jsonlite::unbox(.metric$maximize),
      link_implementation = jsonlite::unbox(get_implementation_url(.config)),
      link_container_image = jsonlite::unbox(get_container_image(.config)),
      component_name = jsonlite::unbox(.config$name),
      references = get_references(.metric),
      additional_info = c(
        get_additional_info(
          .config$info, exclude = c("metrics", "type", "type_info")
        ),
        get_additional_info(
          .metric,
          exclude = c(
            "name",
            "label",
            "summary",
            "description",
            "maximize",
            "min",
            "max",
            "links",
            "references"
          )
        )
      ),
      version = jsonlite::unbox(.config$version)
    )
  })
}) |>
  purrr::list_flatten()

cat("\n>>> Writing output files...\n")
cat("Writing task info to '", par$output, "'...\n", sep = "")
jsonlite::write_json(
  metric_info_json,
  par$output,
  pretty = TRUE,
  null = "null"
)

cat("\n>>> Done!\n")
