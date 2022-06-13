library(RCTD)
library(SingleCellExperiment)
library(Matrix)

# R base for rctd.py

# extract single cell reference data
sce_sc <- sce[, colData(sce)$modality == "sc"]
# get single cell reference counts
sc_counts <- assay(sce_sc, "X")
# get single cell reference labels
sc_cell_types <- factor(colData(sce_sc)$label)
names(sc_cell_types) <- colnames(sce_sc)

# construct reference object (specific for RCTD)
reference <- Reference(sc_counts, sc_cell_types)
# extract spatial data
sce_sp <- sce[, colData(sce)$modality == "sp"]
# get spatial data counts
sp_counts <- assay(sce_sp, "X")
# get spatial data coordinates
sp_coords <- as.data.frame(reducedDim(sce, "spatial"))
colnames(sp_coords) <- c("x", "y")
rownames(sp_coords) <- colnames(sce_sp)
rownames(sp_counts) <- rownames(sce_sp)
colnames(sp_counts) <- colnames(sce_sp)
# create spatial object to use in RCTD
puck <- SpatialRNA(sp_coords, sp_counts)
# create RCTD object from reference and spatialRNA objects
my_rctd <- create.RCTD(puck, reference, max_cores = 1, test_mode = FALSE)
# run analysis and get results
my_rctd <- run.RCTD(my_rctd)
results <- my_rctd@results
cell_type_names <- my_rctd@cell_type_info$info[[2]]
# extract proportions and normalize them (to sum to one)
norm_weights <- sweep(results$weights, 1, rowSums(results$weights), "/")
norm_weights <- as.data.frame(as.matrix(norm_weights))

# proportions estimates needs to be added to colData due
# to bad conversion from reducedDims to obsm (loss of colnames)
colnames(norm_weights) <- paste0("xCT_", cell_type_names)
rownames(norm_weights) <- colnames(sce_sp)
colData(sce_sp) <- cbind(colData(sce_sp), norm_weights)

# return results
sce <- sce_sp

sce
