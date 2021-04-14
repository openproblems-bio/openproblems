## VIASH START
par <- list(
  input = list.files("work", full.names = TRUE, pattern = "*.h5ad"),
  output = "out_bash/modality_alignment/scores.tsv"
)
inp <- par$input[[2]]
## VIASH END

cat("Loading dependencies\n")
library(anndata, warn.conflicts = FALSE)
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(assertthat)

cat("Reading input h5ad files")
scores <- map_df(par$input, function(inp) {
  cat("Reading '", inp, "'\n", sep = "")
  ad <- read_h5ad(inp)

  assert_that("dataset_name" %in% names(ad$uns))
  assert_that("method_name" %in% names(ad$uns))
  assert_that("metric_name" %in% names(ad$uns))
  assert_that("metric_value" %in% names(ad$uns))

  as_tibble(ad$uns[c("dataset_name", "method_name", "metric_name", "metric_value")])
})

write_tsv(scores, par$output)
