import anndata as ad
import h5py
import numpy as np
from dca.api import dca

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
data = np.array(input_file["layers"]["counts"])
input_train = ad.AnnData(data)

print("process data")
# run DCA
dca(input_train, epochs=par["epochs"])


print("Writing data")
# input_train.uns["method_id"] = meta['functionality_name']
# input_train.write_h5ad(par['output'], compression="gzip")
with h5py.File(par['output'], "w") as output_denoised:
    for key in input_file.keys():
        if key == "X":
            continue
        input_file.copy(input_file[key], output_denoised)
    output_denoised["layers"].create_dataset('denoised', data=input_train.X)
    output_denoised["uns"].create_dataset('method_id', data=meta['functionality_name'])

input_file.close()




