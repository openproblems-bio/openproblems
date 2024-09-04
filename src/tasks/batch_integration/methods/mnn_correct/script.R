cat("Loading dependencies\n")
suppressPackageStartupMessages({
  requireNamespace("anndata", quietly = TRUE)
  library(Matrix, warn.conflicts = FALSE)
  requireNamespace("batchelor", quietly = TRUE)
  library(SingleCellExperiment, warn.conflicts = FALSE)
})
## VIASH START
par <- list(
  input = 'resources_test/batch_integration/pancreas/dataset.h5ad',
  output = 'output.h5ad'
)
meta <- list(
  functionality_name = "mnn_correct_feature"
)
## VIASH END

cat("Read input\n")
adata <- anndata::read_h5ad(par$input)

cat("Run mnn\n")
out <- suppressWarnings(batchelor::mnnCorrect(
  t(adata$layers[["normalized"]]),
  batch = adata$obs[["batch"]]
))

cat("Reformat output\n")
layer <- SummarizedExperiment::assay(out, "corrected")
as(t(layer), "sparseMatrix")



cat("Store outputs\n")
output <- anndata::AnnData(
  uns = list(
    dataset_id = adata$uns[["dataset_id"]],
    normalization_id = adata$uns[["normalization_id"]],
    method_id = meta$functionality_name
  ),
  layers = list(
    corrected_counts = as(t(layer), "sparseMatrix")
  ),
  shape = adata$shape
)

cat("Write output to file\n")
zzz <- output$write_h5ad(par$output, compression = "gzip")
