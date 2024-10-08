library(rlang, quietly = TRUE, warn.conflicts = FALSE)
library(purrr, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)

## VIASH START
par <- list(
  "task" = "batch_integration",
  "task_dir" = "src/tasks/batch_integration",
  "output" = "src/tasks/batch_integration/README.md",
  "viash_yaml" = "_viash.yaml",
  "github_url" = "https://github.com/openproblems-bio/openproblems/tree/main/"
)
meta <- list(
  "resources_dir" = "src/common/helper_functions",
  "temp_dir" = "temp/"
)
## VIASH END

if (is.null(par$task) && is.null(par$task_dir)) {
  stop("Either 'task' or 'task_dir' must be provided")
}
if (is.null(par$viash_yaml)) {
  stop("Argument 'viash_yaml' must be provided")
}
if (is.null(par$output)) {
  stop("Argument 'output' must be provided")
}

# import helper function
source(paste0(meta["resources_dir"], "/read_and_merge_yaml.R"))
source(paste0(meta["resources_dir"], "/strip_margin.R"))
source(paste0(meta["resources_dir"], "/read_api_files.R"))

cat("Read task info\n")
task_api <- read_task_api(par[["task_dir"]])

# determine ordering
root <- .task_graph_get_root(task_api)

r_graph <- render_task_graph(task_api, root)

cat("Render API details\n")
order <- names(igraph::bfs(task_api$task_graph, root)$order)
r_details <- map_chr(
  order,
  function(file_name) {
    if (file_name %in% names(task_api$comp_specs)) {
      render_component(task_api$comp_specs[[file_name]])
    } else {
      render_file(task_api$file_specs[[file_name]])
    }
  }
)

cat("Render authors\n")
authors_str <-
  if (nrow(task_api$authors) > 0) {
    paste0(
      "\n## Authors & contributors\n\n",
      task_api$authors %>% knitr::kable() %>% paste(collapse = "\n"),
      "\n"
    )
  } else {
    ""
  }
readme_str <-
  if (is.null(task_api$task_info$readme) || is.na(task_api$task_info$readme)) {
    ""
  } else {
    paste0(
      "\n## README\n\n",
      task_api$task_info$readme,
      "\n"
    )
  }

cat("Generate qmd content\n")
relative_path <- par[["task_dir"]] %>%
  gsub(paste0(dirname(par[["viash_yaml"]]), "/*"), "", .) %>%
  gsub("/*$", "", .)
source_url <- paste0(par[["github_url"]], relative_path)
qmd_content <- strip_margin(glue::glue("
  §---
  §title: \"{task_api$task_info$label}\"
  §format: gfm
  §---
  §
  §<!--
  §This file is automatically generated from the tasks's api/*.yaml files.
  §Do not edit this file directly.
  §-->
  §
  §{task_api$task_info$summary}
  §
  §Path to source: [`{relative_path}`]({source_url})
  §
  §{readme_str}
  §
  §## Motivation
  §
  §{task_api$task_info$motivation}
  §
  §## Description
  §
  §{task_api$task_info$description}
  §{authors_str}
  §## API
  §
  §{r_graph}
  §
  §{paste(r_details, collapse = '\n\n')}
  §
  §"), symbol = "§")

cat("Write README.qmd to file\n")
qmd_file <- tempfile(
  pattern = "README_",
  fileext = ".qmd",
  tmpdir = meta$temp_dir
)

if (!dir.exists(meta$temp_dir)) {
  dir.create(meta$temp_dir, recursive = TRUE)
}
writeLines(qmd_content, qmd_file)

cat("Render README.qmd to README.md\n")
out <- processx::run(
  command = "quarto",
  args = c("render", qmd_file, "--output", "-"),
  echo = TRUE
)

writeLines(out$stdout, par$output)
