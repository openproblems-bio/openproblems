# Dependencies
library(SingleCellExperiment)
library(Seurat)

# Method body
# Default and test parameters
n_pcs <- 50
k_score <- NULL
k_filter <- NULL
if (is_test) {
  n_pcs <- 5
  k_score <- 5
  k_filter <- 20
}

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
