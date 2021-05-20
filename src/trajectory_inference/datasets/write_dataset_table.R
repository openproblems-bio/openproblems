library(httr)
library(tidyverse)

output_dir <- tempfile()
dir.create(output_dir)

# config
deposit_id <- 1443566

print("Retrieveing metadata") 

# retrieve file metadata from zenodo
files <-
  GET(glue::glue("https://zenodo.org/api/records/{deposit_id}")) %>%
  httr::content() %>%
  .$files %>%
  map_df(function(l) {
    as_tibble(t(unlist(l)))
  }) %>%
  filter(grepl("^real/", filename)) %>% 
  mutate(
    name = filename %>% str_replace_all("\\.rds$", "") %>% str_replace_all("/", "_") %>% paste0("zenodo_", deposit_id, "_", .),
    url = paste0("https://github.com/dynverse/dyngen/raw/data_files/", name, ".rds"),
    local_out = paste0(output_dir, "/", name, ".rds")
  )

write.table(files, "src/cellular_dynamics/trajectory_inference/datasets/datasets.tsv")
