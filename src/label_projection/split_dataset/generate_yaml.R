library(tidyverse)
library(anndata)

h5ad_files <- fs::dir_ls("resources/datasets/openproblems_v1", recurse = TRUE, regexp = "\\.h5ad$")

param_list <- map(h5ad_files, function(h5ad_file) {
  ad <- anndata::read_h5ad(h5ad_file, backed = "r")
  if (all(c("batch", "celltype") %in% colnames(ad$obs))) {
    list(
      id = gsub(".*/", "", ad$uns[["dataset_id"]]),
      input = paste0("../../../", h5ad_file)
    )
  } else {
    NULL
  }
})
output <- list(
  param_list = unname(param_list) %>% .[!map_lgl(., is.null)],
  obs_label = "celltype",
  obs_batch = "batch",
  seed = 123L
)

yaml::write_yaml(output, "src/label_projection/split_dataset/params.yaml")
