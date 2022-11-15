import scanpy as sc
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
adata = sc.read_h5ad(par['input'])

print(">> Normalize data")
# libsize and sqrt L1 norm
sqrt_data = scprep.utils.matrix_transform(adata, np.sqrt)
# sqrt_L1, libsize = scprep.normalize.library_size_normalize(sqrt_data, rescale=1, return_library_size=True)
# sqrt_L1 = sqrt_L1.tocsr()

print(sqrt_data)

# print(">> Store output in adata")
# adata.layers["sqrtnorm"] = lognorm
# adata.obs["norm_factor"] = norm["norm_factor"]
# adata.uns["normalization_method"] = meta["functionality_name"].removeprefix("normalize_")

# print(">> Write data")
# adata.write_h5ad(par['output'], compression="gzip")
