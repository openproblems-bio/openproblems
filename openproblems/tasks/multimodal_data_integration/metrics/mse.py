import numpy as np
from scipy import sparse


def mse(adata):
    X = adata.obsm["aligned"]
    Y = adata.uns["mode2"].obsm["aligned"]
    # mean and norm
    Z = np.vstack([X, Y])
    Z -= np.mean(Z, axis=0)
    Z /= np.sqrt(np.sum(Z ** 2))
    # split back out
    X, Y = Z[: X.shape[0]], Z[X.shape[0] :]
    return np.mean(np.sum((X - Y) ** 2))
