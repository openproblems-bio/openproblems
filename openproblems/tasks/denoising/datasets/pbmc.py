import numpy as np
from scipy.sparse import issparse
from scanpy import AnnData
import molecular_cross_validation.util as ut

from ....data.pbmc import load_pbmc
from ....tools.decorators import dataset


def split_data(adata: AnnData, train_frac: float = 0.9, seed: int = 0) -> AnnData:
    """Split data using molecular cross-validation.

    Stores "train" and "test" dataset using the AnnData.obsm property."""

    random_state = np.random.RandomState(seed)

    X = adata.raw.X

    if issparse(X):
        X = np.array(X.todense())
    if np.allclose(X, X.astype(np.int)):
        X = X.astype(np.int)
    else:
        raise TypeError("Molecular cross-validation requires integer count data.")

    X_train, X_test = ut.split_molecules(X, 0.9, 0.0, random_state)
    adata.obsm["train"] = X_train
    adata.obsm["test"] = X_test
    return adata


@dataset("10k PBMCs from a Healthy Donor (10x/v3)")
def pbmc(test=False):
    adata = load_pbmc(test=test)
    adata = split_data(adata)
    return adata
