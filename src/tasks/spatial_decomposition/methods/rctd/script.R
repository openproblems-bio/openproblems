library(anndata)
library(spacexr)
library(Matrix)

## VIASH START
par <- list(
  input_single_cell = "resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/single_cell_ref.h5ad",
  input_spatial_masked = "resources_test/spatial_decomposition/cxg_mouse_pancreas_atlas/spatial_masked.h5ad",
  output = "output.h5ad", 
  fc_cutoff = 0.5, 
  fc_cutoff_reg = 0.75
)
meta <- list(
  functionality_name = "rctd", 
  cpus = 1
)
## VIASH END

cat("Reading input files\n")
input_single_cell <- anndata::read_h5ad(par$input_single_cell)
input_spatial <- anndata::read_h5ad(par$input_spatial)

# set spatial coordinates for the single cell data
coordinates <- matrix(1, dim(input_single_cell)[1], 2)
rownames(coordinates) <- rownames(input_single_cell)
input_single_cell$obsm <- list(coordinates = coordinates)

# remove rare cell types to prevent RCTD error
# celltype_counts <- table(input_single_cell$obs$cell_type)
# input_single_cell <- input_single_cell[input_single_cell$obs$cell_type %in% names(as.table(celltype_counts[celltype_counts > 25]))]

# get single cell reference counts
sc_counts <- t(input_single_cell$layers['counts'])
# get single cell reference labels
sc_cell_types <- factor(input_single_cell$obs$cell_type)
names(sc_cell_types) <- rownames(input_single_cell)
# construct reference object (specific for RCTD)
reference <- Reference(sc_counts, sc_cell_types)

# get spatial data counts
sp_counts <- t(input_spatial$layers['counts'])
# get spatial data coordinates
sp_coords <- as.data.frame(input_spatial$obsm['coordinates'])
colnames(sp_coords) <- c("x", "y")
rownames(sp_coords) <- rownames(input_spatial)
# create spatial object to use in RCTD
puck <- SpatialRNA(sp_coords, sp_counts)

# create RCTD object from reference and spatialRNA objects
if (!is.null(meta$cpus)) {
max_cores <- meta$cpus
} else {
max_cores <- 1
}
rctd <- create.RCTD(
  puck,
  reference,
  max_cores = max_cores,
  fc_cutoff = par$fc_cutoff,
  fc_cutoff_reg = par$fc_cutoff_reg,
  test_mode = FALSE,
  UMI_min_sigma = 100,
  CELL_MIN_INSTANCE = 1
)

# run analysis and get results
rctd <- run.RCTD(rctd)
results <- rctd@results
cell_type_names <- rctd@cell_type_info$info[[2]]

# extract proportions and normalize them (to sum to one)
norm_weights <- sweep(results$weights, 1, rowSums(results$weights), "/")
norm_weights <- as.matrix(norm_weights)
coordinates <- as.matrix(sp_coords)

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
    coordinates = coordinates,
    proportions_pred = norm_weights
  ),
  layers = list(
    counts = input_spatial$layers['counts']
  )
)
output$write_h5ad(par[["output"]], compression = "gzip")
