import anndata as ad
import scprep
import numpy as np

## VIASH START
par = {
    'input': "output_train.h5ad",
    'output': "output_norm.h5ad"
}
meta = {
    "functionality_name": "normalize_sqrt_L1"
}
## VIASH END

print(">> Load data")
adata = ad.read_h5ad(par['input'])

print(">> Normalize data")
# libsize and sqrt L1 norm
sqrt_data = scprep.utils.matrix_transform(adata.X, np.sqrt)
sqrt_L1, libsize = scprep.normalize.library_size_normalize(sqrt_data, rescale=1, return_library_size=True)
sqrt_L1 = sqrt_L1.tocsr()

print(">> Store output in adata")
adata.layers["sqrtnorm"] = sqrt_L1
adata.uns["normalization_method"] = meta["functionality_name"].removeprefix("normalize_")

print(adata.to_df(layer="sqrtnorm"))

print(">> Write data")
adata.write_h5ad(par['output'], compression="gzip")
