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

  if ("normalization_method" %in% names(ad$uns)) {
    uns_names <- c("dataset_id", "normalization_method", "method_id", "metric_id", "metric_value")
  } else {
    uns_names <- c("dataset_id", "method_id", "metric_id", "metric_value")
  }

  for (uns_name in uns_names) {
    assert_that(
      uns_name %in% names(ad$uns),
      msg = paste0("File ", inp, " must contain `uns['", uns_name, "']`")
    )
  }

  as_tibble(ad$uns[uns_names])
})

write_tsv(scores, par$output)
