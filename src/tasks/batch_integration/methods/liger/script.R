cat(">> Load dependencies\n")
requireNamespace("anndata", quietly = TRUE)
requireNamespace("rliger", quietly = TRUE)

## VIASH START
par <- list(
  input = "resources_test/batch_integration/pancreas/dataset.h5ad",
  output = "output.h5ad"
)
meta <- list(
  functionality_name = "liger"
)
## VIASH END

cat("Read input\n")
adata <- anndata::read_h5ad(par$input)

anndataToLiger <- function(adata) {
  # fetch batch names
  batch <- adata$obs$batch
  batch_names <- as.character(unique(batch))

  # restructure data
  raw_data <- lapply(batch_names, function(batch_name) {
    Matrix::t(adata$layers[["counts"]][batch == batch_name, , drop = FALSE])
  })
  names(raw_data) <- batch_names

  rliger::createLiger(raw.data = raw_data, remove.missing = FALSE)
}

addNormalizedDataToLiger <- function(adata, lobj) {
  norm_data <- lapply(names(lobj@raw.data), function(name) {
    norm <- adata$layers[["normalized"]]
    # subset
    norm <- norm[
      colnames(lobj@raw.data[[name]]),
      rownames(lobj@raw.data[[name]]),
      drop = FALSE
    ]
    # transpose
    norm <- Matrix::t(norm)

    # turn into dgcMatrix
    as(as(norm, "denseMatrix"), "CsparseMatrix")
  })
  names(norm_data) <- names(lobj@raw.data)

  lobj@norm.data <- norm_data

  lobj
}

cat(">> Create Liger Data object\n")
lobj <- anndataToLiger(adata)

cat(">> Normalize data\n")
lobj <- addNormalizedDataToLiger(adata, lobj)

# could also use the rliger normalization instead
# lobj <- rliger::normalize(lobj)

cat(">> Select genes\n")
# lobj <- rliger::selectGenes(lobj)
# overwrite gene selection to include all genes
lobj@var.genes <- adata$var_names

cat(">> Perform scaling\n")
lobj <- rliger::scaleNotCenter(lobj, remove.missing = FALSE)

cat(">> Joint Matrix Factorization\n")
lobj <- rliger::optimizeALS(lobj, k = 20)

cat(">> Quantile normalization\n")
lobj <- rliger::quantile_norm(lobj)

cat(">> Store dimred in adata\n")
adata$obsm[["X_emb"]] <- lobj@H.norm[rownames(adata), , drop = FALSE]
adata$uns[["method_id"]] <- meta$functionality_name

cat(">> Write AnnData to disk\n")
zzz <- adata$write_h5ad(par$output, compression = "gzip")
