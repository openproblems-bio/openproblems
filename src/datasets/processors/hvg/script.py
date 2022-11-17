
import scanpy as sc

### VIASH START
par = {
  'input': 'work/ca/0751ff85df6f9478cb7bda5a705cad/zebrafish.sqrt_cpm.pca.output.h5ad',
  'layer_input': 'normalized',
  'output': 'dataset.h5ad',
  'var_hvg': 'hvg',
  'var_hvg_score': 'hvg_score',
  'num_features': 100
}
### VIASH END

print(">> Load data")
adata = sc.read(par['input'])

print(">> Look for layer")
layer = adata.X if not par['layer_input'] else adata.layers[par['layer_input']]

print(">> Run PCA")
out = sc.pp.highly_variable_genes(
  adata,
  layer=par["layer_input"],
  n_top_genes=par["num_features"],
  flavor='cell_ranger',
  inplace=False
)

print(">> Storing output")
adata.var[par["var_hvg"]] = out['highly_variable'].values
adata.var[par["var_hvg_score"]] = out['dispersions'].values

print(">> Writing data")
adata.write_h5ad(par['output'])

