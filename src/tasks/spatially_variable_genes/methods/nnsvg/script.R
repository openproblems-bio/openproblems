suppressMessages(library(SpatialExperiment))
suppressMessages(library(scran))
suppressMessages(library(nnSVG))
suppressMessages(library(anndata))
suppressMessages(library(dplyr))

# VIASH START
par = list(
    'input_data' = 'resources_test/spatially_variable_genes/mouse_brain_coronal/dataset.h5ad',
    'output' = 'output.h5ad'
)
meta = list(
    'functionality_name' = 'nnSVG',
    'cpus' = 4
)

# VIASH END

# load data
cat('Load data')
adata <- read_h5ad(par$input_data)
counts <- t(as.matrix(adata$layers[['counts']]))
    
colnames(counts) <- adata$obs_names
rownames(counts) <- adata$var_names
    
loc <- as.data.frame(adata$obsm[['spatial']])

row_data = adata$var
row_data$gene_id = rownames(row_data)
row_data$feature_type = "Gene Expression"

colnames(loc) <- c("x", "y")
rownames(loc) <- colnames(counts)

spe <- SpatialExperiment(
    assays = list(counts = counts),
    rowData = row_data,
    colData = loc, 
    spatialCoordsNames = c("x", "y"))

# calculate logcounts (log-transformed normalized counts) using scran package
# using library size factors
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)

# run nnSVG
if (!is.null(meta$cpus)) {
    n_cpus <- meta$cpus
} else {
    n_cpus <- 1
}

cat('Run nnSVG')
spe <- nnSVG(spe, n_threads=n_cpus)

df <- as.data.frame(rowData(spe)) %>%
    subset(select = c('feature_id', 'LR_stat'))

colnames(df) <- c('feature_id', 'pred_spatial_var_score')
rownames(df) <- NULL

# save output
cat("Write output AnnData to file\n")
output = anndata::AnnData(
    shape = adata$shape, 
    var=df,
    uns=list('dataset_id' = adata$uns[['dataset_id']],
             'method_id' =  meta[['functionality_name']]))

anndata::write_h5ad(anndata = output, filename = par$output)
