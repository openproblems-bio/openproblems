import anndata as ad
import h5py
import numpy as np
from dca.api import dca
from scipy import sparse

## VIASH START
par = {
    'input_train': 'output/denoising/pancreas_split_data_output_train.h5ad',
    'output': 'output_dca.h5ad',
    'epochs': 300,
}
meta = {
    'functionality_name': 'dca',
}
## VIASH END

print("load input data")
# input_train = ad.read_h5ad(par['input_train'])
input_file = h5py.File(par["input_train"], 'r')

# convert to csr matrix
group = input_file.get("layers/counts")
matrix = sparse.csr_matrix((group["data"], group["indices"], group["indptr"]))

# Convert to Anndata
input_train = ad.AnnData(matrix, dtype=matrix.dtype)

print("process data")
# run DCA
dca(input_train, epochs=par["epochs"])

# Convert to csr matrix
output_train = sparse.csr_matrix(input_train.X)

print("Writing data")


with h5py.File(par['output'], "w") as output_denoised:
    for key in input_file.keys():
        input_file.copy(input_file[key], output_denoised)
    denoised = output_denoised['layers'].create_group('denoised')
    denoised.create_dataset('data', data=output_train.data)
    denoised.create_dataset('indices', data=output_train.indices)
    denoised.create_dataset('indptr', data=output_train.indptr)
    for key, value in dict(input_file['layers/counts'].attrs).items():
        output_denoised['layers/denoised'].attrs[key] = value
    output_denoised["uns"].create_dataset('method_id', data=meta['functionality_name'])

input_file.close()




