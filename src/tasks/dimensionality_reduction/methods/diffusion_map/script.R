requireNamespace("anndata", quietly = TRUE)
requireNamespace("diffusionMap", quietly = TRUE)

## VIASH START
par <- list(
  input = "resources_test/dimensionality_reduction/pancreas/dataset.h5ad",
  output = "output.h5ad",
  n_dim = 3
)
## VIASH END

cat("Reading input files\n")
input <- anndata::read_h5ad(par$input)

cat("Running destiny diffusion map\n")
# create SummarizedExperiment object
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(
    logcounts = t(as.matrix(input$layers[["normalized"]]))
  )
)
dm <- destiny::DiffusionMap(sce)
X_emb <- destiny::eigenvectors(dm)[, seq_len(par$n_dim)]

cat("Write output AnnData to file\n")
output <- anndata::AnnData(
  uns = list(
    dataset_id = input$uns[["dataset_id"]],
    normalization_id = input$uns[["normalization_id"]],
    method_id = meta$functionality_name
  ),
  obsm = list(
    X_emb = X_emb
  ),
  shape = input$shape
)
output$write_h5ad(par$output, compression = "gzip")
