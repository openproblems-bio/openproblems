import anndata as ad
from numba import njit
from typing import Tuple
import numpy as np
from sklearn.metrics import pairwise_distances
from scipy.sparse import issparse


_K = 30

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


@njit(cache=True, fastmath=True)
def _coranking_matrix(R1: np.ndarray, R2: np.ndarray) -> np.ndarray:  # pragma: no cover
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


@njit(cache=True, fastmath=True)
def _continuity(Q: np.ndarray, m: int) -> np.ndarray:  # pragma: no cover

    C = np.zeros(m - 1)  # continuity

    for k in range(m - 1):
        Qs = Q[:k, k:]
        # a row vector of weights. weight = rank error = actual_rank - k
        W = np.arange(Qs.shape[1]).reshape(1, -1)
        # 1 - normalized hard-k-extrusions. upper-right region
        C[k] = 1 - np.sum(Qs * W) / ((k + 1) * m * (m - 1 - k))

    return C


@njit(cache=True, fastmath=True)
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


## VIASH START
par = {
    'input_reduced': 'resources_test/dimensionality_reduction/pancreas/reduced.h5ad',
    'input_test': 'resources_test/dimensionality_reduction/pancreas/test.h5ad',
    'output': 'score.h5ad',
}
meta = {
    'functionality_name': 'nn_ranking',
}
## VIASH END

print("Load data")
input_reduced = ad.read_h5ad(par['input_reduced'])
input_test = ad.read_h5ad(par['input_test'])

# Select 1000 most variable genes
idx = input_test.var['hvg_score'].to_numpy().argsort()[-1000:]
input_test = input_test[:, idx]

# Compute pairwise distances
if issparse(input_test):
    Dx = pairwise_distances(input_test.layers['normalized'].A)
else:
    Dx = pairwise_distances(input_test.layers['normalized'])

De = pairwise_distances(input_reduced.obsm["X_emb"])
Rx, Re = _ranking_matrix(Dx), _ranking_matrix(De)
Q = _coranking_matrix(Rx, Re)

T, C, QNN, AUC, LCMC, _kmax, Qlocal, Qglobal = _metrics(Q)

print("Store metric value")
input_reduced.uns['metric_ids'] =  {meta['functionality_name']: ['continuity', 'co-KNN size', 'co-KNN AUC', 'local continuity meta criterion', 'local property', 'global property']}
if np.any(np.isnan(input_reduced.obsm["X_emb"])):
    input_reduced.uns['metric_values'] = [0.0, 0.0, 0.0, 0.5, -np.inf, -np.inf, -np.inf]
else:
    input_reduced.uns['metric_values'] = [C[_K], QNN[_K], AUC, LCMC[_K], Qlocal, Qglobal]


print("Copy data to new AnnData object")
output = ad.AnnData(
    uns={}
)
output.uns['normalization_id'] = input_reduced.uns['normalization_id']
output.uns['method_id'] = input_reduced.uns['method_id']
output.uns['dataset_id'] = input_reduced.uns['dataset_id']
output.uns['metric_ids'] =  input_reduced.uns['metric_ids']
output.uns['metric_values'] = input_reduced.uns['metric_values']

print("Write data to file")
output.write_h5ad(par['output'], compression="gzip")