cat(">> Loading dependencies\n")
library(anndata, warn.conflicts = FALSE)
library(ALRA, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input_train = "resources_test/denoising/pancreas/train.h5ad",
  norm = "log",
  output = "output.h5ad"
)
meta <- list(
  functionality_name = "alra"
)
## VIASH END

cat(">> Load input data\n")
input_train <- read_h5ad(par$input_train, backed = "r")

cat(">> Set normalization method\n")
if (par$norm == "sqrt") {
  norm_fn <- sqrt
  denorm_fn <- function(x) x^2
} else if (par$norm == "log") {
  norm_fn <- log1p
  denorm_fn <- expm1
} else {
  stop("Unknown normalization method: ", par$norm)
}

cat(">> Normalize data\n")
data <- as.matrix(input_train$layers[["counts"]])
totalPerCell <- rowSums(data)
data <- sweep(data, 1, totalPerCell, "/")
data <- norm_fn(data)

cat(">> Run ALRA\n")
data <- alra(data)$A_norm_rank_k_cor_sc
data <- denorm_fn(data)
data <- sweep(data, 1, totalPerCell, "*")

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
