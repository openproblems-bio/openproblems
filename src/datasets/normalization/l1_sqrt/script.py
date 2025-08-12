import anndata as ad
import scprep
import numpy as np

## VIASH START
par = {
    'input': "output_train.h5ad",
    'output': "output_norm.h5ad"
}
meta = {
    'name': "l1_sqrt"
}
## VIASH END

print("Load data", flush=True)
adata = ad.read_h5ad(par['input'])

print("Normalize data", flush=True)
# libsize and sqrt L1 norm
sqrt_data = scprep.utils.matrix_transform(adata.layers['counts'], np.sqrt)
l1_sqrt, libsize = scprep.normalize.library_size_normalize(sqrt_data, rescale=1, return_library_size=True)
l1_sqrt = l1_sqrt.tocsr()

print("Store output in adata", flush=True)
adata.layers[par["layer_output"]] = l1_sqrt
adata.uns["normalization_id"] = par["normalization_id"] or meta['name']

print("Write data", flush=True)
adata.write_h5ad(par['output'], compression="gzip")
