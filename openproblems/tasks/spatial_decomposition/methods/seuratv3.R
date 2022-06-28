library(Seurat)

# R base for seuratv3.py
sce_sc <- as.Seurat(sce_sc, counts = "X", data = NULL)
sce_sp <- as.Seurat(sce_sp, counts = "X", data = NULL)

print(sce_sp)

# Normalize and do dimred for spatial data
sce_sp <- SCTransform(sce_sp, assay = "originalexp", verbose = FALSE)

sce_sp <- RunPCA(sce_sp, assay = "SCT", verbose = FALSE)
sce_sp <- RunUMAP(sce_sp, reduction = "pca", dims = 1:30)

# Normalize and do dimred for single cell data
sce_sc <- SCTransform(
  sce_sc,
  assay = "originalexp", ncells = min(3000, nrow(sce_sc)), verbose = FALSE
)
sce_sc <- RunPCA(sce_sc, verbose = FALSE)
sce_sc <- RunUMAP(sce_sc, dims = 1:30)

# find anchors (MNN's to compute adjustmen vectors)
anchors <- FindTransferAnchors(
  reference = sce_sc,
  query = sce_sp,
  normalization.method = "SCT"
)

# transfer labels from single cell data to spatial
predictions_assay <- TransferData(
  anchorset = anchors,
  refdata = as.factor(as.character(sce_sc@meta.data$label)),
  prediction.assay = TRUE,
  weight.reduction = sce_sp[["pca"]],
  dims = 1:30
)

# format data and return results
predictions <- GetAssayData(predictions_assay, slot = "data")
predictions <- predictions[!(rownames(predictions) == "max"), ]
rownames(predictions) <- paste0("xCT_", rownames(predictions))
colnames(predictions) <- colnames(sce_sp)
predictions <- as.data.frame(t(predictions))
sce_sp@meta.data <- cbind(sce_sp@meta.data, predictions)
sce_sp <- as.SingleCellExperiment(sce_sp)

sce_sp
