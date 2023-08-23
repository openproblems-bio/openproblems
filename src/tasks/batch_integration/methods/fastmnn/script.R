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
  output = 'output.h5ad'
)
meta <- list(
  functionality_name = "mnn_correct_feature"
)
## VIASH END

cat("Read input\n")
adata <- anndata::read_h5ad(par$input)

# TODO: pass output of 'multiBatchNorm' to fastMNN

cat("Run mnn\n")
out <- suppressWarnings(batchelor::fastMNN(
  t(adata$layers[["normalized"]]),
  batch = adata$obs[["batch"]]
))

cat("Reformat output\n")
obsm <- SingleCellExperiment::reducedDim(out, "corrected")
adata$obsm[["X_emb"]] <- obsm
# return_type == "feature" is currently not working in fastMNN

# # reusing the same script for mnn_correct and mnn_correct_feature
# return_type <- gsub("mnn_correct_", "", meta[["functionality_name"]])

# if (return_type == "feature") {
#   layer <- SummarizedExperiment::assay(out, "corrected")
#   adata$layers[["corrected_counts"]] <- as(t(layer), "sparseMatrix")
# } else if (return_type == "embedding") {
#   obsm <- SingleCellExperiment::reducedDim(out, "corrected")
#   adata$obsm[["X_emb"]] <- obsm
# }

cat("Store outputs\n")
adata$uns[["method_id"]] <- meta$functionality_name
zzz <- adata$write_h5ad(par$output, compression = "gzip")
