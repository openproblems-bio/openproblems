import sys
import anndata as ad
import numpy as np

## VIASH START
par = {
    'input': "resources_test/common/pancreas/dataset.h5ad",
    'output_train': "train.h5ad",
    'output_test': "test.h5ad",
    'train_frac': 0.9,
    'seed': 0
}
meta = {
    "functionality_name": "process_dataset",
    "resources_dir": "src/tasks/denoising/process_dataset"
}
## VIASH END

# add helper scripts to path
sys.path.append(meta["resources_dir"])
from helper import split_molecules

# set random state
random_state = np.random.RandomState(par['seed'])

print(">> Load Data", flush=True)
adata = ad.read_h5ad(par["input"])

# remove all layers except for counts
for key in list(adata.layers.keys()):
    if key != "counts":
        del adata.layers[key]

# round counts and convert to int
counts = np.array(adata.layers["counts"]).round().astype(int)

print(">> process and split data", flush=True)
train_data, test_data = split_molecules(
    counts.data, par["train_frac"], 0.0, random_state
)

X_train = counts.copy()
X_test = counts.copy()
X_train.data = train_data
X_test.data = test_data
X_train.eliminate_zeros()
X_test.eliminate_zeros()

# copy adata to train_set, test_set
output_train = ad.AnnData(
    layers={"counts": X_train},
    obs=adata.obs[[]],
    var=adata.var[[]],
    uns={"dataset_id": adata.uns["dataset_id"]}
)
test_uns_keys = ["dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description", "dataset_organism"]
output_test = ad.AnnData(
    layers={"counts": X_test},
    obs=adata.obs[[]],
    var=adata.var[[]],
    uns={key: adata.uns[key] for key in test_uns_keys}
)

# add additional information for the train set
output_test.uns["train_sum"] = X_train.sum()

# Remove no cells that do not have enough reads
is_missing = np.array(X_train.sum(axis=0) == 0)

output_train = output_train[:, ~is_missing.flatten()]
output_test = output_test[:, ~is_missing.flatten()]

print(">> Write to file", flush=True)
output_train.write_h5ad(par["output_train"])
output_test.write_h5ad(par["output_test"])
