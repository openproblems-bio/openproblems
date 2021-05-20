## VIASH START
par <- list(
  output =  "datasets",
  download_link = ""
)
## VIASH END

library(httr)
library(tidyverse)

output_dir <- tempfile()
dir.create(output_dir)

tmp <- tempfile()
on.exit(file.remove(tmp))

download.file(download_link, tmp, quiet = TRUE)

ds <- read_rds(tmp)
ad <- to_h5ad(ds)
ad$uns["dataset_id"] <- ds$id

write_h5ad(ad, paste0(par$output + "/" + ad$uns["dataset_id"]))

