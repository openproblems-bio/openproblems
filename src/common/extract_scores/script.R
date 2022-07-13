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

  for (uns_name in par$column_names) {
    assert_that(
      uns_name %in% names(ad$uns),
      msg = paste0("File ", inp, " must contain `uns['", uns_name, "']`")
    )
  }

  as_tibble(ad$uns[uns_names])
})

write_tsv(scores, par$output)
