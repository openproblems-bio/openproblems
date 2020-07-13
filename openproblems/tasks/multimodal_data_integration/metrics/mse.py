import numpy as np
from scipy import sparse
import scprep


def _square(X):
    if sparse.issparse(X):
        X.data = X.data ** 2
        return X
    else:
        return scprep.utils.toarray(X) ** 2


def mse(adata):
    X = adata.obsm["aligned"]
    Y = adata.obsm["mode2_aligned"]
    # mean and norm
    vstack = sparse.vstack if sparse.issparse(X) and sparse.issparse(Y) else np.vstack
    Z = vstack([X, Y])
    Z -= np.mean(Z, axis=0)
    Z /= np.sqrt(np.sum(_square(Z)))
    # split back out
    X, Y = Z[: X.shape[0]], Z[X.shape[0] :]
    return np.mean(np.sum(_square(X - Y)))
