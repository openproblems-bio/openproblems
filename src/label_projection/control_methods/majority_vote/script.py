import anndata as ad

## VIASH START
par = {
    'input_train': 'resources_test/label_projection/pancreas/train.h5ad',
    'input_test': 'resources_test/label_projection/pancreas/test.h5ad',
    'output': 'output.h5ad'
}
meta = {
    'functionality_name': 'foo'
}
## VIASH END

print("Load data")
input_train = ad.read_h5ad(par['input_train'])
input_test = ad.read_h5ad(par['input_test'])

print("Compute majority vote")
majority = input_train.obs.label.value_counts().index[0]

print("Create prediction object")
input_test.obs["label_pred"] = majority

print("Write output to file")
input_test.uns["method_id"] = meta["functionality_name"]
input_test.write_h5ad(par["output"], compression="gzip")
