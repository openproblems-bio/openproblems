import scanpy as sc
import numpy as np
import scipy.sparse
import molecular_cross_validation.util

## VIASH START
par = {
    'input': "/home/kai/Documents/openroblems/openproblems-v2/resources_test/common/pancreas/dataset.h5ad",
    'output_train': "output/nextflow/pancreas_split_data_output_train.h5d",
    'output_test': "output/nextflow/pancreas_split_data_output_test.h5d",
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
adata = sc.read_h5ad(par["input"])


# remove all layers except for counts
for key in list(adata.layers.keys()):
    if key != "counts":
        del adata.layers[key]

adata.X = adata.layers["counts"]
X = np.array(adata.X)

# for test purposes
X = X.round()

print(">> process and split data")
if scipy.sparse.issparse(X):
    X = np.array(X.todense())
if np.allclose(X, X.astype(int)):
    X = X.astype(int)
else:
    raise TypeError("Molecular cross-validation requires integer count data.")

X_train, X_test = molecular_cross_validation.util.split_molecules(
    X, par["train_frac"], 0.0, random_state
)

# Remove no cells that do not have enough reads
is_missing = X_train.sum(axis=0) == 0
X_train, X_test = X_train[:, ~is_missing], X_test[:, ~is_missing]

#   copy adata to train_set, test_set

output_train = adata[:, ~is_missing].copy()
del output_train.layers["counts"]
output_train.layers["counts"] = scipy.sparse.csr_matrix(X_train).astype(float)

output_test = adata[:, ~is_missing].copy()
del output_test.layers["counts"]
output_test.layers["counts"] = scipy.sparse.csr_matrix(X_test).astype(float)


# NOTE: remove cells which have not enough reads -> is it possible to remove same cells/genes in test set as in train set a above ?
#   * first, remove genes with 0 expression
#   * then, remove cells with 1 count or less

# def filter_genes_cells(adata):
#     """Remove empty cells and genes."""
#     sc.pp.filter_genes(adata, min_cells=1)
#     sc.pp.filter_cells(adata, min_counts=1)
#     return adata


# filter_genes_cells(output_train)

print(">> Writing")

output_train.write_h5ad(par["output_train"])
output_test.write_h5ad(par["output_test"])
