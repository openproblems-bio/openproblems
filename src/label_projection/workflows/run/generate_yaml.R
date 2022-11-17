library(tidyverse)
library(anndata)

h5ad_files <- fs::dir_ls("resources/label_projection/openproblems_v1", recurse = TRUE, regexp = "\\.h5ad$")

regex <- ".*/([^\\.]*)\\.([^\\.]*)\\.([^\\.]*)\\.([^\\.]*)\\.h5ad"

df <- tibble(
  path = as.character(h5ad_files),
  dataset_id = gsub(regex, "\\1", path),
  normalization_id = gsub(regex, "\\2", path),
  comp = gsub(regex, "\\3", path),
  arg_name = gsub(regex, "\\4", path)
) %>%
  spread(arg_name, path)

param_list <- pmap(df, function(dataset_id, normalization_id, comp, output_solution, output_test, output_train) {
  list(
    id = paste0(dataset_id, ".", normalization_id),
    input_train = output_train,
    input_test = output_test,
    input_solution = output_solution,
    dataset_id = dataset_id,
    normalization_id = normalization_id
  )
})

output <- list(
  param_list = param_list
  # obs_label = "celltype",
  # obs_batch = "batch",
  # seed = 123L
)

yaml::write_yaml(output, "src/label_projection/workflows/run/params.yaml")
