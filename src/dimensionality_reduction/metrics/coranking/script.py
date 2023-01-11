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

import anndata as ad
from numba import njit
from typing import Tuple
import numpy as np
from sklearn.metrics import pairwise_distances
from scipy.sparse import issparse


_K = 30

## VIASH START
par = {
    "input_reduced": "resources_test/dimensionality_reduction/pancreas/reduced.h5ad",
    "input_test": "resources_test/dimensionality_reduction/pancreas/test.h5ad",
    "output": "score.h5ad",
}
## VIASH END

print("Load data", flush=True)
input_reduced = ad.read_h5ad(par["input_reduced"])
input_test = ad.read_h5ad(par["input_test"])

X_emb = input_reduced.obsm["X_emb"]
high_dim = input_test.layers["normalized"]

def _ranking_matrix(pdist_mat: np.ndarray) -> np.ndarray:
    """The pairwise distance matrix is ranked using argsort."""
    assert pdist_mat.shape[0] == pdist_mat.shape[1]
    return np.argsort(np.argsort(pdist_mat))

def _coranking_matrix(rmat1: np.ndarray, rmat2: np.ndarray) -> np.ndarray:
    """Compute the coranking matrix from two ranking matrices."""
    assert rmat1.shape == rmat2.shape
    m = rmat1.shape[0]
    corank_mat = np.zeros((m - 1, m - 1), dtype=np.int32)
    for i in range(m):
        for j in range(m):
            if i != j:
                k = rmat1[i, j] - 1
                l = rmat2[i, j] - 1
                corank_mat[k, l] += 1
    return corank_mat

def _metrics(
    Q: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float, np.ndarray, int, float, float]:
    Q = Q[1:, 1:]
    m = len(Q)

    T = _trustworthiness(Q, m)
    C = _continuity(Q, m)
    QNN = _qnn(Q, m)
    LCMC = _lcmc(QNN, m)
    kmax = _kmax(LCMC)
    Qlocal = _q_local(QNN, kmax)
    Qglobal = _q_global(QNN, kmax, m)
    AUC = _qnn_auc(QNN)

    return T, C, QNN, AUC, LCMC, kmax, Qlocal, Qglobal

def _fit(
    original: np.ndarray, embedding: np.ndarray
) -> Tuple[float, float, float, float, float, float, float]:
    from sklearn.metrics import pairwise_distances

    if np.any(np.isnan(E)):
        return 0.0, 0.0, 0.0, 0.5, -np.inf, -np.inf, -np.inf

    pdist_original = pairwise_distances(original)
    pdist_embedding = pairwise_distances(embedding)
    rmat_original = _ranking_matrix(pdist_original)
    rmat_embedding = _ranking_matrix(pdist_embedding)
    corank = _coranking_matrix(rmat_original, rmat_embedding)

    T, C, QNN, AUC, LCMC, _kmax, Qlocal, Qglobal = _metrics(corank)

    return T[_K], C[_K], QNN[_K], AUC, LCMC[_K], Qlocal, Qglobal


print("Store metric value", flush=True)
input_reduced.uns['metric_ids'] =  {meta['functionality_name']: ['continuity', 'co-KNN size', 'co-KNN AUC', 'local continuity meta criterion', 'local property', 'global property']}
if np.any(np.isnan(input_reduced.obsm["X_emb"])):
    input_reduced.uns['metric_values'] = [0.0, 0.0, 0.0, 0.5, -np.inf, -np.inf, -np.inf]
else:
    input_reduced.uns['metric_values'] = [C[_K], QNN[_K], AUC, LCMC[_K], Qlocal, Qglobal]


print("Copy data to new AnnData object", flush=True)
output = ad.AnnData(
    uns={key: input_reduced.uns[key] for key in ["dataset_id", "normalization_id", "method_id", 'metric_ids', 'metric_values']}
)

print("Write data to file", flush=True)
output.write_h5ad(par['output'], compression="gzip")




















def _trustworthiness(Q: np.ndarray, m: int) -> np.ndarray:  # pragma: no cover

    T = np.zeros(m - 1)  # trustworthiness

    for k in range(m - 1):
        Qs = Q[k:, :k]
        # a column vector of weights. weight = rank error = actual_rank - k
        W = np.arange(Qs.shape[0]).reshape(-1, 1)
        # 1 - normalized hard-k-intrusions. lower-left region.
        # weighted by rank error (rank - k)
        T[k] = 1 - np.sum(Qs * W) / ((k + 1) * m * (m - 1 - k))

    return T


def _continuity(Q: np.ndarray, m: int) -> np.ndarray:  # pragma: no cover

    C = np.zeros(m - 1)  # continuity

    for k in range(m - 1):
        Qs = Q[:k, k:]
        # a row vector of weights. weight = rank error = actual_rank - k
        W = np.arange(Qs.shape[1]).reshape(1, -1)
        # 1 - normalized hard-k-extrusions. upper-right region
        C[k] = 1 - np.sum(Qs * W) / ((k + 1) * m * (m - 1 - k))

    return C


def _qnn(Q: np.ndarray, m: int) -> np.ndarray:  # pragma: no cover

    QNN = np.zeros(m)  # Co-k-nearest neighbor size

    for k in range(m):
        # Q[0,0] is always m. 0-th nearest neighbor is always the point itself.
        # Exclude Q[0,0]
        QNN[k] = np.sum(Q[: k + 1, : k + 1]) / ((k + 1) * m)

    return QNN


def _lcmc(QNN: np.ndarray, m: int) -> np.ndarray:
    LCMC = QNN - (np.arange(m) + 1) / (m - 1)
    return LCMC


def _kmax(LCMC: np.ndarray) -> int:
    kmax = np.argmax(LCMC)
    return kmax  # type: ignore


def _q_local(QNN: np.ndarray, kmax: int) -> float:
    Qlocal = np.sum(QNN[: kmax + 1]) / (kmax + 1)
    return Qlocal


def _q_global(QNN: np.ndarray, kmax: int, m: int) -> float:
    # skip the last. The last is (m-1)-nearest neighbor, including all samples.
    Qglobal = np.sum(QNN[kmax:-1]) / (m - kmax - 1)
    return Qglobal


def _qnn_auc(QNN: np.ndarray) -> float:
    AUC = np.mean(QNN)
    return AUC  # type: ignore
