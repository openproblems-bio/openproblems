library(tidyverse)
library(anndata)

h5ad_files <- fs::dir_ls("resources/label_projection/openproblems_v1", recurse = TRUE, regexp = "\\.h5ad$")

regex <- ".*/([^\\.]*)\\.([^\\.]*)\\.([^\\.]*)\\.h5ad"

df <- tibble(
  path = as.character(h5ad_files),
  id = gsub(regex, "\\1", path),
  comp = gsub(regex, "\\2", path),
  arg_name = gsub(regex, "\\3", path)
) %>%
  spread(arg_name, path)

param_list <- pmap(df, function(id, comp, output_solution, output_test, output_train) {
  list(
    id = id,
    input_train = output_train,
    input_test = output_test,
    input_solution = output_solution
  )
})

output <- list(
  param_list = param_list
  # obs_label = "celltype",
  # obs_batch = "batch",
  # seed = 123L
)

yaml::write_yaml(output, "src/label_projection/workflows/run/params.yaml")
