import anndata as ad
import numpy as np
import scprep
from magic import MAGIC
import scipy


## VIASH START
par = {
    'input_train': 'output_train.h5ad',
    'output': 'output_magic.h5ad',
    'solver': 'exact',
    'norm': 'sqrt',
    'decay': 1,
    't': 3,
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
    input_train.layers['counts'], rescale=1, return_library_size=True
)

X = scprep.utils.matrix_transform(X, norm_fn)
Y = MAGIC(solver=par['solver'], verbose=False, decay=par['decay'], t=par['t']).fit_transform(
    X, genes="all_genes"
)

Y = scprep.utils.matrix_transform(Y, denorm_fn)
Y = scprep.utils.matrix_vector_elementwise_multiply(Y, libsize, axis=0)

output_denoised = input_train.copy()
output_denoised.uns["method_id"] = meta["functionality_name"]
output_denoised.layers["denoised"] = scipy.sparse.csr_matrix(Y)

print("Writing Data")
output_denoised.write_h5ad(par['output'],compression="gzip")

