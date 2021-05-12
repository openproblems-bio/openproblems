## VIASH START
par <- list(
  output =  "datasets"
)
## VIASH END


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

print("Downloading datasets & converting to AnnData & writing to file") 

# iterate over rows
pbapply::pblapply(
  seq_len(nrow(files)),
  cl = 8,
  function(i) {
    tmp <- tempfile()
    on.exit(file.remove(tmp))
    
    # download file
    download.file(files$links.download[[i]], tmp, quiet = TRUE)
    
    # read counts
    ds <- read_rds(tmp)
    ad <- to_h5ad(ds)
    ad$uns["dataset_id"] = ds$id
    
    write_h5ad(ad, paste0(par$output + "/" + ad$uns["dataset_id"]))
  }
)

