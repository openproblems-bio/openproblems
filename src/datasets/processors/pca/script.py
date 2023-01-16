
import scanpy as sc

### VIASH START
par = {
  'input': 'resources_test/common/pancreas/dataset.h5ad',
  'layer_input': 'log_cpm',
  'output': 'dataset.h5ad',
  'obsm_embedding': 'X_pca',
  'varm_loadings': 'pca_loadings',
  'uns_variance': 'pca_variance',
  'num_components': 25
}
### VIASH END

print(">> Load data")
adata = sc.read(par['input'])

print(">> Look for layer")
layer = adata.X if not par['layer_input'] else adata.layers[par['layer_input']]

print(">> Run PCA")
X_pca, loadings, variance, variance_ratio = sc.tl.pca(
    layer, 
    n_comps=par["num_components"], 
    return_info=True
)

print(">> Storing output")
adata.obsm[par["obsm_embedding"]] = X_pca
adata.varm[par["varm_loadings"]] = loadings.T
adata.uns[par["uns_variance"]] = {
    "variance": variance, 
    "variance_ratio": variance_ratio
}

print(">> Writing data")
adata.write_h5ad(par['output'])

