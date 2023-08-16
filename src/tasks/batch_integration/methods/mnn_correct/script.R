cat("Loading dependencies\n")
suppressPackageStartupMessages({
  requireNamespace("anndata", quietly = TRUE)
  library(Matrix, warn.conflicts = FALSE)
  requireNamespace("batchelor", quietly = TRUE)
  library(SingleCellExperiment, warn.conflicts = FALSE)
})
## VIASH START
par <- list(
  input = 'resources_test/batch_integration/pancreas/unintegrated.h5ad',
  output = 'output.h5ad',
  hvg = FALSE
)
meta <- list(
  functionality_name = "mnn_correct_feature"
)
## VIASH END

cat("Read input\n")
adata <- anndata::read_h5ad(par$input)

# don't subset when return_type is not "feature"
if ("hvg" %in% names(par) && par$hvg) {
  cat("Select HVGs\n")
  adata <- adata[, adata$var[["hvg"]]]
}

cat("Run mnn\n")
out <- suppressWarnings(batchelor::mnnCorrect(
  t(adata$layers[["normalized"]]),
  batch = adata$obs[["batch"]]
))

cat("Reformat output\n")
layer <- SummarizedExperiment::assay(out, "corrected")
adata$layers[["corrected_counts"]] <- as(t(layer), "sparseMatrix")

cat("Store outputs\n")
adata$uns[["method_id"]] <- meta$functionality_name
zzz <- adata$write_h5ad(par$output, compression = "gzip")
