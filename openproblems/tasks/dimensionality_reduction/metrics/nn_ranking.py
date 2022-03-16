"""
This file is uses slightly modified code from pyDRMetrics [1]_, see:

    - https://doi.org/10.1016/j.heliyon.2021.e06199 - the article.
    - https://data.mendeley.com/datasets/jbjd5fmggh/1 - the supplementary files.

The following changes have been made:

    - :mod:`numba` JIT for performance reasons
    - use broadcasting instead of a 3rd loop in :func:`_ranking_matrix`

[1] Zhang, Yinsheng (2021),
  “Source code, sample data, and case study report for pyDRMetrics”,
  Mendeley Data, V1, doi: 10.17632/jbjd5fmggh.1
"""

from ....tasks.dimensionality_reduction.methods.preprocessing import preprocess_scanpy
from ....tools.decorators import metric
from anndata import AnnData
from numba import njit
from scipy.sparse import issparse
from sklearn.metrics import pairwise_distances
from typing import Tuple

import numpy as np

__original_author__ = "Yinsheng Zhang"
__original_author_email__ = "zhangys@illinois.edu"
__license__ = "CC BY 4.0"
__license_link__ = (
    "https://data.mendeley.com/datasets/"
    "jbjd5fmggh/1/files/da1bca42-c4da-4376-9177-bd2d9a308108"
)


_K = 30


@njit(cache=True, fastmath=True)
def _ranking_matrix(D: np.ndarray) -> np.ndarray:
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


@njit(cache=True, fastmath=True)
def _coranking_matrix(R1: np.ndarray, R2: np.ndarray) -> np.ndarray:
    assert R1.shape == R2.shape
    Q = np.zeros(R1.shape, dtype=np.int32)
    m = len(Q)
    for i in range(m):
        for j in range(m):
            k = int(R1[i, j])
            l = int(R2[i, j])  # noqa: E741
            Q[k, l] += 1

    return Q


@njit(cache=True, fastmath=True)
def _metrics(
    Q: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, np.ndarray, int, float, float]:
    Q = Q[1:, 1:]
    m = len(Q)

    T = np.zeros(m - 1)  # trustworthiness
    C = np.zeros(m - 1)  # continuity
    QNN = np.zeros(m)  # Co-k-nearest neighbor size
    LCMC = np.zeros(m)  # Local Continuity Meta Criterion

    for k in range(m - 1):
        Qs = Q[k:, :k]
        # a column vector of weights. weight = rank error = actual_rank - k
        W = np.arange(Qs.shape[0]).reshape(-1, 1)
        # 1 - normalized hard-k-intrusions. lower-left region.
        # weighted by rank error (rank - k)
        T[k] = 1 - np.sum(Qs * W) / (k + 1) / m / (m - 1 - k)
        Qs = Q[:k, k:]
        # a row vector of weights. weight = rank error = actual_rank - k
        W = np.arange(Qs.shape[1]).reshape(1, -1)
        # 1 - normalized hard-k-extrusions. upper-right region
        C[k] = 1 - np.sum(Qs * W) / (k + 1) / m / (m - 1 - k)

    for k in range(m):
        # Q[0,0] is always m. 0-th nearest neighbor is always the point itself.
        # Exclude Q[0,0]
        QNN[k] = np.sum(Q[: k + 1, : k + 1]) / ((k + 1) * m)
        LCMC[k] = QNN[k] - (k + 1) / (m - 1)

    kmax = np.argmax(LCMC)
    Qlocal = np.sum(QNN[: kmax + 1]) / (kmax + 1)
    # skip the last. The last is (m-1)-nearest neighbor, including all samples.
    Qglobal = np.sum(QNN[kmax:-1]) / (m - kmax - 1)
    AUC = np.mean(QNN)

    return T, C, QNN, AUC, LCMC, kmax, Qlocal, Qglobal


def _high_dim(adata: AnnData) -> np.ndarray:
    preprocess_scanpy(adata)
    high_dim = adata.X
    return high_dim.A if issparse(high_dim) else high_dim


def _fit(
    X: np.ndarray, E: np.ndarray
) -> Tuple[float, float, float, float, float, float, float]:
    if np.any(np.isnan(E)):
        return 0.0, 0.0, 0.0, 0.5, -np.inf, -np.inf, -np.inf

    Dx = pairwise_distances(X)
    De = pairwise_distances(E)
    Rx, Re = _ranking_matrix(Dx), _ranking_matrix(De)
    Q = _coranking_matrix(Rx, Re)

    T, C, QNN, AUC, LCMC, _kmax, Qlocal, Qglobal = _metrics(Q)

    return T[_K], C[_K], QNN[_K], AUC, LCMC[_K], Qlocal, Qglobal


@metric("continuity", maximize=True)
def continuity(adata: AnnData) -> float:
    _, C, _, *_ = _fit(_high_dim(adata), adata.obsm["X_emb"])
    return float(np.clip(C, 0.0, 1.0))  # in [0, 1]


@metric("co-KNN size", maximize=True)
def qnn(adata: AnnData) -> float:
    _, _, QNN, *_ = _fit(_high_dim(adata), adata.obsm["X_emb"])
    # normalized in the code to [0, 1]
    return float(np.clip(QNN, 0.0, 1.0))


@metric("co-KNN AUC", maximize=True)
def qnn_auc(adata: AnnData) -> float:
    _, _, _, AUC, *_ = _fit(_high_dim(adata), adata.obsm["X_emb"])
    return float(np.clip(AUC, 0.5, 1.0))  # in [0.5, 1]


@metric("local continuity meta criterion", maximize=True)
def lcmc(adata: AnnData) -> float:
    *_, LCMC, _, _ = _fit(_high_dim(adata), adata.obsm["X_emb"])
    return LCMC


@metric("local property", maximize=True)
def qlocal(adata: AnnData) -> float:
    # according to authors, this is usually preferred to
    # qglobal, because human are more sensitive to nearer neighbors
    *_, Qlocal, _ = _fit(_high_dim(adata), adata.obsm["X_emb"])
    return Qlocal


@metric("global property", maximize=True)
def qglobal(adata: AnnData) -> float:
    *_, Qglobal = _fit(_high_dim(adata), adata.obsm["X_emb"])
    return Qglobal
