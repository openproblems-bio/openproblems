import anndata as ad
import numpy as np

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

print("Load data", flush=True)
input_train = ad.read_h5ad(par['input_train'])
input_test = ad.read_h5ad(par['input_test'])

print("Compute label distribution", flush=True)
label_distribution = input_train.obs.label.value_counts()
label_distribution = label_distribution / label_distribution.sum()

print("Create prediction object", flush=True)
input_test.obs["label_pred"] = np.random.choice(
    label_distribution.index,
    size=input_test.n_obs,
    replace=True,
    p=label_distribution
)

print("Write output to file", flush=True)
input_test.uns["method_id"] = meta["functionality_name"]
input_test.write_h5ad(par["output"], compression="gzip")
