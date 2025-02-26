import scanpy as sc

# Load your processed data
adata = sc.read_h5ad('output.h5ad')

# Basic information
print(adata)

# Check metadata
print(adata.uns['dataset_id'])
print(adata.uns['dataset_name'])

# Look at cell annotations
print(adata.obs.head())

# Check available layers
print(list(adata.layers.keys()))