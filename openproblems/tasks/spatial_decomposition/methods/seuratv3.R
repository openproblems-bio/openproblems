#' Seuratv3 TransferData
#' @param sce_sc SingleCellExperiment single-cell data
#' @param sce_sp SingleCellExperiment spatial data
#' @param n_pcs int Number of principal components
#' @param sctransform_n_cells int Number of cells sampled to build NB regression

options(error = rlang::entrace)

library(Seurat)
library(future)

# 8GB max size up from default of 500MB
options(future.globals.maxSize = 8 * 1024^3)
plan(multicore, workers = availableCores())

args <- readRDS("/tmp/openproblems_seurat_args.rds")
sce_sc <- args$sce_sc
sce_sp <- args$sce_sp
n_pcs <- args$n_pcs

# R base for seuratv3.py
sce_sc <- as.Seurat(sce_sc, counts = "X", data = NULL)
sce_sp <- as.Seurat(sce_sp, counts = "X", data = NULL)

# Normalize and do dimred for spatial data
sce_sp <- SCTransform(
  sce_sp,
  assay = "originalexp",
  ncells = min(sctransform_n_cells, nrow(sce_sp)),
  verbose = TRUE
)

sce_sp <- RunPCA(sce_sp, assay = "SCT", verbose = FALSE, n_pcs = n_pcs)

# Normalize and do dimred for single cell data
sce_sc <- SCTransform(
  sce_sc,
  assay = "originalexp",
  ncells = min(sctransform_n_cells, nrow(sce_sc)),
  verbose = TRUE
)
sce_sc <- RunPCA(sce_sc, verbose = FALSE, n_pcs = n_pcs)

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
  dims = 1:n_pcs
)

# format data and return results
predictions <- GetAssayData(predictions_assay, slot = "data")
predictions <- predictions[!(rownames(predictions) == "max"), ]
rownames(predictions) <- paste0("xCT_", rownames(predictions))
colnames(predictions) <- colnames(sce_sp)
predictions <- as.data.frame(t(predictions))
sce_sp@meta.data <- cbind(sce_sp@meta.data, predictions)
sce_sp <- as.SingleCellExperiment(sce_sp)

saveRDS(sce_sp, "/tmp/openproblems_seurat_sce_sp_out.rds")
