library(rlang, quietly = TRUE, warn.conflicts = FALSE)
library(purrr, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)

## VIASH START
par <- list(
  "task" = "batch_integration",
  "output" = "src/tasks/batch_integration/README.md",
  "viash_yaml" = "_viash.yaml"
)
meta <- list(
  "resources_dir" = "src/common/helper_functions",
  "temp_dir" = "temp/"
)
## VIASH END

# import helper function
source(paste0(meta["resources_dir"], "/read_and_merge_yaml.R"))
source(paste0(meta["resources_dir"], "/strip_margin.R"))
source(paste0(meta["resources_dir"], "/read_api_files.R"))

cat("Read task info\n")
task_dir <- paste0(dirname(par[["viash_yaml"]]), "/src/tasks/", par[["task"]]) %>%
  gsub("^\\./", "", .)
task_api <- read_task_api(task_dir)

r_graph <- render_task_graph(task_api)

# todo: fix hard coded node
order <- names(igraph::bfs(task_api$task_graph, "file_common_dataset")$order)

cat("Render API details\n")
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

cat("Generate qmd content\n")
task_dir_short <- gsub(".*openproblems-v2/", "", task_dir)
qmd_content <- strip_margin(glue::glue("
  §---
  §title: \"{task_api$task_info$label}\"
  §format: gfm
  §---
  §
  §{task_api$task_info$summary}
  §
  §Path: [`{task_dir_short}`](https://github.com/openproblems-bio/openproblems-v2/tree/main/{task_dir_short})
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
md_content <- system(
  paste0("quarto render ", qmd_file, " --output -"),
  ignore.stderr = TRUE,
  intern = TRUE
)

writeLines(md_content, par$output)