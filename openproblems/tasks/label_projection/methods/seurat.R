# Dependencies
library(SingleCellExperiment)
library(Seurat)

# Method body
# Default and test parameters
n.pcs <- 50
k.score <- NULL
k.filter <- NULL
if (is_test) {
  n.pcs <- 5
  k.score <- 5
  k.filter <- 20
}

data <- as.Seurat(sce, counts = "X", data = NULL)
data <- SCTransform(data, assay = "originalexp", verbose = FALSE)

reference <- data[, data$is_train]
query <- data[, !data$is_train]

reference <- RunPCA(reference, assay = "SCT", verbose = FALSE, n_pcs = n.pcs)
anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = "SCT",
  dims = 1:n.pcs,
  k.score = k.score,
  k.filter = k.filter,
  verbose = FALSE
)
query <- TransferData(
  anchorset = anchors,
  reference = reference,
  query = query,
  k.weight = k.score,
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
