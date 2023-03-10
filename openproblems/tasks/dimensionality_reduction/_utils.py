from numba import njit

import numpy as np


@njit(cache=True, fastmath=True)
def _ranking_matrix(D: np.ndarray) -> np.ndarray:  # pragma: no cover
    assert D.shape[0] == D.shape[1]
    R = np.zeros(D.shape)
    m = len(R)
    ks = np.arange(m)

    for i in range(m):
        for j in range(m):
            R[i, j] = np.sum(
                (D[i, :] < D[i, j]) | ((ks < j) & (np.abs(D[i, :] - D[i, j]) <= 1e-12))
            )

    return R


def ranking_matrix(X):
    from sklearn.metrics import pairwise_distances

    D = pairwise_distances(X)
    R = _ranking_matrix(D)
    return R
