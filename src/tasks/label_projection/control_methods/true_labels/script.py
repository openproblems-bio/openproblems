import anndata as ad

## VIASH START
par = {
    'input_train': 'resources_test/label_projection/pancreas/train.h5ad',
    'input_test': 'resources_test/label_projection/pancreas/test.h5ad',
    'input_solution': 'resources_test/label_projection/pancreas/test.h5ad',
    'output': 'output.h5ad'
}
meta = {
    'functionality_name': 'foo'
}
## VIASH END

print("Load data", flush=True)
# input_train = ad.read_h5ad(par['input_train'])
input_test = ad.read_h5ad(par['input_test'])
input_solution = ad.read_h5ad(par['input_solution'])

print("Create prediction object", flush=True)
input_test.obs["label_pred"] = input_solution.obs["label"]

print("Write output to file", flush=True)
input_test.uns["method_id"] = meta["functionality_name"]
input_test.write_h5ad(par["output"], compression="gzip")
