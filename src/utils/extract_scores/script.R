## VIASH START
par <- list(
  input = list.files("out_bash/modality_alignment/metrics/", full.names = TRUE),
  output = "out_bash/modality_alignment/scores.tsv"
)
inp <- par$input[[1]]
## VIASH END

cat("Loading dependencies\n")
library(anndata, warn.conflicts = FALSE)
options(tidyverse.quiet = TRUE)
library(tidyverse)

cat("Reading input h5ad files")
scores <- map_df(par$input, function(inp) {
  ad <- read_h5ad(inp)
  as_tibble(ad$uns[c("dataset_name", "method_name", "metric_name", "metric_value")])
})

write_tsv(scores, par$output)
