#' seurat.R
#' Runs Seurat TransferData with SCT
#' @param sce SingleCellExperiment
#' @param n_pcs int Number of PCA components
#' @param k_filter int How many neighbors (k) to use when filtering anchors.
#' Set to NA to turn off filtering.
#' @param k_score int How many neighbors (k) to use when scoring anchors

# Dependencies
library(SingleCellExperiment)
library(Seurat)
library(Matrix)

# Convert RsparseMatrix
if (is(assay(sce, "X"), "RsparseMatrix") {
  assay(sce, "X") <- as(assay(sce, "X"), "CsparseMatrix")
}

# Method body

data <- as.Seurat(sce, counts = "X", data = NULL)
data <- SCTransform(data, assay = "originalexp", verbose = FALSE)

reference <- data[, data$is_train]
query <- data[, !data$is_train]

reference <- RunPCA(reference, assay = "SCT", verbose = FALSE, n_pcs = n_pcs)
anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = "SCT",
  dims = 1:n_pcs,
  k.score = k_score,
  k.filter = k_filter,
  verbose = FALSE
)
query <- TransferData(
  anchorset = anchors,
  reference = reference,
  query = query,
  k.weight = k_score,
  refdata = list(
    labels = "labels"
  ),
  verbose = FALSE
)
query$labels_pred <- query$predicted.labels
reference$labels_pred <- reference$labels
sce <- as.SingleCellExperiment(merge(reference, query))

# Return
sce
