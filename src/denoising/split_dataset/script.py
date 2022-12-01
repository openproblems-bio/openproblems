import anndata as ad
import numpy as np
import scipy.sparse
import molecular_cross_validation.util

## VIASH START
par = {
    'input': "resources_test/common/pancreas/dataset.h5ad",
    'output_train': "train.h5ad",
    'output_test': "test.h5ad",
    'train_frac': 0.9,
    'seed': 0
}
meta = {
    "functionality_name": "split_data"
}
## VIASH END
"""Split data using molecular cross-validation.

Stores "train" and "test" dataset in separate ad files.
"""

random_state = np.random.RandomState(par['seed'])

print(">> Load Data")
adata = ad.read_h5ad(par["input"])


# remove all layers except for counts
for key in list(adata.layers.keys()):
    if key != "counts":
        del adata.layers[key]

counts_rounded = np.array(adata.layers["counts"]).round()

counts = counts_rounded.astype(int)

print(">> process and split data")
train_data, test_data = molecular_cross_validation.util.split_molecules(
    counts.data, par["train_frac"], 0.0, random_state
)

X_train = counts.copy()
X_test = counts.copy()
X_train.data = train_data
X_test.data = test_data

#   copy adata to train_set, test_set
output_train = ad.AnnData(
    layers={"counts": X_train.astype(float)},
    obs=adata.obs[[]],
    var=adata.var[[]],
    uns={"dataset_id": adata.uns["dataset_id"]}
)
output_test = ad.AnnData(
    layers={"counts": X_test.astype(float)},
    obs=adata.obs[[]],
    var=adata.var[[]],
    uns={"dataset_id": adata.uns["dataset_id"]}
)

# Remove no cells that do not have enough reads
is_missing = np.array(X_train.sum(axis=0) == 0)

output_train = output_train[:, ~is_missing.flatten()]
output_test = output_test[:, ~is_missing.flatten()]

# output_test = adata[:, ~is_missing].copy()
# del output_test.layers["counts"]
# output_test.layers["counts"] = scipy.sparse.csr_matrix(X_test).astype(float)


print(">> Writing")
output_train.write_h5ad(par["output_train"])
output_test.write_h5ad(par["output_test"])
