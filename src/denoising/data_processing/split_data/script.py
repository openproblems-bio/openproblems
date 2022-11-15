import scanpy as sc
import numpy as np
import scipy.sparse
import molecular_cross_validation.util

## VIASH START
par = {
    'input': "resources_test/common/pancreas/dataset.h5ad",
    'output_train': "output_train.h5ad",
    'output_test': "output_test.h5ad",
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


X = adata.layers["counts"]

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
    X, 0.9, 0.0, random_state
)



# copy adata to train_set, test_set

output_train = adata
output_train.layers["counts"] = scipy.sparse.csr_matrix(X_train).astype(float)

output_test = adata
output_test.layers["counts"] = scipy.sparse.csr_matrix(X_test).astype(float)

# TODO: remove zero entries -> uncertain how this is done. below code gives error that matrix size is different
# is_missing = output_train.layers["counts"].sum(axis=0) == 0
# output_train.layers["counts"], output_test.layers["counts"] = output_train.layers["counts"][:, ~is_missing], output_test.layers["counts"][:, ~is_missing]

print(">> Writing")

output_train.write_h5ad(par["output_train"])
output_test.write_h5ad(par["output_test"])
