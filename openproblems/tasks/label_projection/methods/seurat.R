# Dependencies
library(SingleCellExperiment)
library(Seurat)

# Method body
data <- as.Seurat(sce, counts="X")
data <- SCTransform(data, verbose = FALSE)
reference <- data[, data$is_train]
query <- data[, !data$is_train]
anchors <- FindTransferAnchors(
  reference = reference,
  query = query,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)
query <- TransferData(
  anchorset = anchors, 
  reference = reference,
  query = query,
  refdata = list(
    labels = "labels",
  )
)
query$labels_pred <- query$predicted.labels
reference$labels_pred <- reference$labels
sce <- merge(reference, query)

# Return
sce