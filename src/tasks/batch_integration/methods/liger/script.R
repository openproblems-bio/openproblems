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

  rliger::createLiger(rawData = raw_data, removeMissing = FALSE)
}

addNormalizedDataToLiger <- function(adata, lobj) {
  norm_data <- lapply(names(rliger::rawData(lobj)), function(name) {
    norm <- adata$layers[["normalized"]]

    # subset
    col_names <- colnames(rliger::rawData(lobj)[[name]])
    row_names <- rownames(rliger::rawData(lobj)[[name]])
    prefix <- paste0(name, "_")
    col_names <- sub(prefix, "", col_names)

    norm <- norm[
      col_names,
      row_names,
      drop = FALSE
    ]

    # add prefix
    rownames(norm) <- paste0(prefix, rownames(norm))

    # transpose
    norm <- Matrix::t(norm)

    # turn into dgcMatrix
    as(as(norm, "denseMatrix"), "CsparseMatrix")
  })
  names(norm_data) <- names(rliger::rawData(lobj))

  for (name in names(rliger::rawData(lobj))) {
    lobj@datasets[[name]]@normData <- norm_data[[name]]
  }

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
lobj@varFeatures <- adata$var_names

cat(">> Perform scaling\n")
lobj <- rliger::scaleNotCenter(lobj, removeMissing = FALSE)

cat(">> Joint Matrix Factorization\n")
lobj <- rliger::runIntegration(lobj, k = 20)

cat(">> Quantile normalization\n")
lobj <- rliger::quantileNorm(lobj)

cat(">> Store output\n")
# remove dataset names from rownames
for (name in names(rliger::rawData(lobj))) {
  rownames(lobj@H.norm) <- sub(paste0(name, "_"), "", rownames(lobj@H.norm))
}

output <- anndata::AnnData(
  uns = list(
    dataset_id = adata$uns[["dataset_id"]],
    normalization_id = adata$uns[["normalization_id"]],
    method_id = meta$functionality_name
  ),
  obsm = list(
    X_emb = lobj@H.norm[rownames(adata), , drop = FALSE]
  ),
  shape = adata$shape
)

cat(">> Write AnnData to file\n")
zzz <- output$write_h5ad(par$output, compression = "gzip")
