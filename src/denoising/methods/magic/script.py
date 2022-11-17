import anndata as ad
import numpy as np
import scprep
from magic import MAGIC


## VIASH START
par = {
    'input_train': 'output_train.h5ad',
    'output': 'output_magic.h5ad',
    'layer_input': 'counts',
    'solver': 'exact',
    'norm': 'sqrt'
}
meta = {
    'functionality_name': 'foo',
}
## VIASH END


print("load data")
input_train = ad.read_h5ad(par['input_train'])

normtype = par['norm']

if normtype == "sqrt":
    norm_fn = np.sqrt
    denorm_fn = np.square
elif normtype == "log":
    norm_fn = np.log1p
    denorm_fn = np.expm1

print("processing data")

X, libsize = scprep.normalize.library_size_normalize(
    input_train.layers[par['layer_input']], rescale=1, return_library_size=True
)

X = scprep.utils.matrix_transform(X, norm_fn)
Y = MAGIC(solver=par['solver'], verbose=False).fit_transform(
    X, genes="all_genes"
)

Y = scprep.utils.matrix_transform(Y, denorm_fn)
Y = scprep.utils.matrix_vector_elementwise_multiply(Y, libsize, axis=0)

output_denoised = input_train.copy()
output_denoised.uns["method_id"] = meta["functionality_name"]
output_denoised.layers["denoised"] = Y

print("Writing Data")
output_denoised.write_h5ad(par['output'],compression="gzip")

