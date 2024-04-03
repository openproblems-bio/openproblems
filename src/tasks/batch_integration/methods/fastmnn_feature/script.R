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
# reusing the same script for fastmnn_embed and fastmnn_feature
return_type <- gsub("fastmnn_", "", meta[["functionality_name"]])

output <- anndata::AnnData(
  shape = adata$shape,
  uns = list(
    dataset_id = adata$uns[["dataset_id"]],
    normalization_id = adata$uns[["normalization_id"]],
    method_id = meta$functionality_name
  )
)

if (return_type == "feature") {
  layer <- as(SummarizedExperiment::assay(out, "reconstructed"), "sparseMatrix")
  output$layers[["corrected_counts"]] <- t(layer)
} else if (return_type == "embedding") {
  obsm <- SingleCellExperiment::reducedDim(out, "corrected")
  output$obsm[["X_emb"]] <- obsm
}

cat("Write output to file\n")
zzz <- output$write_h5ad(par$output, compression = "gzip")
