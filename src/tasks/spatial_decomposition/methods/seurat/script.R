library(anndata)
library(Seurat)

## VIASH START
par <- list(
  input_single_cell = "resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/single_cell_ref.h5ad",
  input_spatial_masked = "resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/spatial_masked.h5ad",
  output = "output.h5ad", 
  n_pcs = 30,
  sctransform_n_cells = 500
)
meta <- list(
  functionality_name = "seurat"
)
## VIASH END

cat("Reading input files\n")
input_single_cell <- anndata::read_h5ad(par$input_single_cell)
input_spatial <- anndata::read_h5ad(par$input_spatial)

cat(">> Converting AnnData to Seurat\n")
anndataToSeurat <- function(adata, assay) {
  obj <- SeuratObject::CreateSeuratObject(counts = as(Matrix::t(adata$layers[["counts"]]), "CsparseMatrix"), assay = assay)
  obj <- SeuratObject::AddMetaData(object = obj, metadata = adata$obs)
  obj
}

seurat_sc <- anndataToSeurat(input_single_cell, "RNA")
seurat_sp <- anndataToSeurat(input_spatial, "spatial")

cat(">> Generate predictions\n")

# Normalize and do dimred for spatial data
seurat_sp <- SCTransform(
  seurat_sp,
  assay = "spatial",
  ncells = min(par$sctransform_n_cells, nrow(seurat_sp)),
  verbose = TRUE,
  conserve.memory = TRUE
)

seurat_sp <- RunPCA(seurat_sp, assay = "SCT", verbose = FALSE, n_pcs = par$n_pcs)

# Normalize and do dimred for single cell data
seurat_sc <- SCTransform(
  seurat_sc,
  assay = "RNA",
  ncells = min(par$sctransform_n_cells, nrow(seurat_sc)),
  verbose = TRUE,
  conserve.memory = TRUE
)

seurat_sc <- RunPCA(seurat_sc, verbose = FALSE, n_pcs = par$n_pcs)

# find anchors (MNN's to compute adjustmen vectors)
anchors <- FindTransferAnchors(
  reference = seurat_sc,
  query = seurat_sp,
  normalization.method = "SCT"
)

# transfer labels from single cell data to spatial
predictions_assay <- TransferData(
  anchorset = anchors,
  refdata = as.factor(as.character(seurat_sc@meta.data$cell_type)),
  prediction.assay = TRUE,
  weight.reduction = seurat_sp[["pca"]],
  dims = 1:par$n_pcs
)

# format data and return results
predictions <- LayerData(predictions_assay, layer = "data")
predictions <- predictions[!(rownames(predictions) == "max"), ]
predictions <- t(predictions)

sp_coords <- as.data.frame(input_spatial$obsm['coordinates'])
colnames(sp_coords) <- c("x", "y")
rownames(sp_coords) <- rownames(input_spatial)
sp_coords <- as.matrix(sp_coords)

cat("Write output AnnData to file\n")
output <- anndata::AnnData(
  shape = input_spatial$shape, 
  obs = input_spatial$obs,
  var = input_spatial$var,
  uns = list(
    cell_type_names = input_spatial$uns['cell_type_names'],
    dataset_id = input_spatial$uns[["dataset_id"]],
    method_id = meta[["functionality_name"]]
  ),
  obsm = list(
    coordinates = sp_coords,
    proportions_pred = predictions
  ),
  layers = list(
    counts = input_spatial$layers['counts']
  )
)
output$write_h5ad(par[["output"]], compression = "gzip")
