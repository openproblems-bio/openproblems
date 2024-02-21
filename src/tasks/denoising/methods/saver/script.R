cat(">> Loading dependencies\n")
library(anndata, warn.conflicts = FALSE)
library(SAVER, warn.conflicts = FALSE)
library(Matrix, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input_train = "resources_test/denoising/pancreas/train.h5ad",
  norm = "log",
  output = "output.h5ad"
)
meta <- list(
  functionality_name = "saver",
  ncpus = 30
)
## VIASH END

cat(">> Load input data\n")
input_train <- read_h5ad(par$input_train, backed = "r")

cat(">> Normalize data\n")
data <- as(t(input_train$layers[["counts"]]), "CsparseMatrix")

cat(">> Run SAVER\n")
data <- t(saver(data, ncores = meta$ncpus, estimates.only = TRUE))

cat(">> Store output\n")
output <- AnnData(
  layers = list(denoised = data),
  obs = input_train$obs[, c(), drop = FALSE],
  var = input_train$var[, c(), drop = FALSE],
  uns = list(
    dataset_id = input_train$uns[["dataset_id"]],
    method_id = meta$functionality_name
  )
)

cat(">> Write output to file\n")
output$write_h5ad(par$output, compression = "gzip")
