library(httr)
library(tidyverse)

output_dir <- tempfile()
dir.create(output_dir)

# config
deposit_id <- 1443566

cat("Retrieving metadata\n")

# retrieve file metadata from zenodo
files <-
  GET(glue::glue("https://zenodo.org/api/records/{deposit_id}")) %>%
  httr::content() %>%
  .$files %>%
  map_df(function(l) {
    as_tibble(t(unlist(l)))
  }) %>%
  transmute(
    id = filename %>% str_replace_all("\\.rds$", "") %>% str_replace_all("/", "_") %>% paste0("zenodo_", deposit_id, "_", .),
    url = links.download,
    checksum,
    filesize
  )

# subsample dataset
files <- files %>% slice(1:10)

# TODO rename links.download as header to links_download
write_tsv(files, "src/trajectory_inference/datasets/download_datasets/datasets.tsv")
