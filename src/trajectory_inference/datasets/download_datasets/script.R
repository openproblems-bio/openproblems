## VIASH START
par <- list(
  id = "ti_dataset",
  output =  "dataset",
  input1 = "https://zenodo.org/api/files/8b17ae8e-2fd8-4ab6-9b3b-a1def87cdf34/real/silver/placenta-trophoblast-differentiation-invasive_mca.rds"
)
## VIASH END

print("> Loading dependencies")

options(tidyverse.quiet = TRUE) # make sure tidyverse is quiet
library(tidyverse)
requireNamespace("dynio", quietly = TRUE)
requireNamespace("anndata", quietly = TRUE)

output_dir <- tempfile()
dir.create(output_dir)

cat("> Checking input parameter\n")

input_path <-
  if (grepl("^https?://", par$input)) {
    cat("> Downloading file from remote\n")
    # Check if link or local file
    tmp <- tempfile()
    on.exit(file.remove(tmp))
    utils::download.file(par$input, tmp, quiet = TRUE)
    tmp
  } else {
    par$input
  }

cat("> Reading file\n")
ds <- read_rds(input_path)

cat("> Converting RDS to h5ad\n")
ad <- dynio::to_h5ad(ds)

cat("> Writing to h5ad\n")
ad$write_h5ad(paste0(par$output))
