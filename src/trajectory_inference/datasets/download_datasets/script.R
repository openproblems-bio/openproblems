## VIASH START
par <- list(
  id = "ti_dataset",
  output =  "dataset",
  input1 = "https://zenodo.org/api/files/8b17ae8e-2fd8-4ab6-9b3b-a1def87cdf34/real/silver/placenta-trophoblast-differentiation-invasive_mca.rds"
)
## VIASH END

print("1")

library(httr)
library(tidyverse)
library(dynio)
library(anndata)

output_dir <- tempfile()
dir.create(output_dir)

print("2")

if(par$input1 %>% startsWith("http://") || par$input1 %>% startsWith("https://")){
  # Check if link or local file
  tmp <- tempfile()
  on.exit(file.remove(tmp))
  
  download.file(par$input1, tmp, quiet = TRUE)
  
} else {
  print("3")
  tmp <- par$input1
}

ds <- read_rds(tmp)
ad <- to_h5ad(ds)
print("Here")
ad$write_h5ad(paste0(par$output))
