
import scanpy as sc

### VIASH START
par = {
  'input': 'work/ca/0751ff85df6f9478cb7bda5a705cad/zebrafish.sqrt_cpm.pca.output.h5ad',
  'layer_input': 'normalized',
  'output': 'dataset.h5ad',
  'key_added': 'knn',
  'n_neighbors': 15
}
### VIASH END

print(">> Load data")
adata = sc.read(par['input'])

print(">> Look for layer")
adata.X = adata.layers[par['layer_input']]

print(">> Run kNN")
sc.pp.neighbors(
    adata,
    use_rep='X_pca',
    key_added=par['key_added'],
    n_neighbors=par['num_neighbors']
)

del adata.X

print(">> Writing data")
adata.write_h5ad(par['output'])

