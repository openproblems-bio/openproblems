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
    X = scprep.utils.toarray(adata.obsm["aligned"])
    Y = scprep.utils.toarray(adata.obsm["mode2_aligned"])

    # mean and norm
    Z = np.vstack([X, Y])
    Z -= np.mean(Z, axis=0)
    Z /= np.sqrt(np.sum(_square(Z)))

    # split back out
    X, Y = Z[: X.shape[0]], Z[X.shape[0] :]
    return np.mean(np.sum(_square(X - Y)))
