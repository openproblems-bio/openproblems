import numpy as np
from scipy import sparse
import scprep


def mse(adata):
    X = scprep.utils.to_array_or_spmatrix(adata.obsm["aligned"])
    Y = scprep.utils.to_array_or_spmatrix(adata.uns["mode2"].obsm["aligned"])
    # mean and norm
    vstack = sparse.vstack if sparse.issparse(X) and sparse.issparse(Y) else np.vstack
    Z = vstack([X, Y])
    Z -= np.mean(Z, axis=0)
    Z /= np.sqrt(np.sum(Z ** 2))
    # split back out
    X, Y = Z[: X.shape[0]], Z[X.shape[0] :]
    return np.mean(np.sum((X - Y) ** 2))
