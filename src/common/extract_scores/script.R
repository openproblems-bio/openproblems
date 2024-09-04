cat("Loading dependencies\n")
library(anndata, warn.conflicts = FALSE)
options(tidyverse.quiet = TRUE)
library(tidyverse)
library(assertthat)

## VIASH START
par <- list(
  input = "resources_test/label_projection/pancreas/knn_accuracy.h5ad",
  output = "scores.tsv"
)
inp <- par$input[[1]]
## VIASH END

cat("Reading input h5ad files\n")
scores <- map_df(par$input, function(inp) {
  cat("Reading '", inp, "'\n", sep = "")
  ad <- read_h5ad(inp)

  for (uns_name in par$column_names) {
    assert_that(
      uns_name %in% names(ad$uns),
      msg = paste0("File ", inp, " must contain `uns['", uns_name, "']`")
    )
  }

  data.frame(ad$uns[par$column_names])
})

write_tsv(scores, par$output)
