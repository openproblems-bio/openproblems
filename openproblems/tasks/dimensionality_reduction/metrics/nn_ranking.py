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

from ....tools.decorators import metric
from ..utils import ranking_matrix
from anndata import AnnData
from numba import njit
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


def _fit(adata: AnnData) -> Tuple[float, float, float, float, float, float, float]:
    Rx = adata.obsm["X_ranking"]
    E = adata.obsm["X_emb"]

    if np.any(np.isnan(E)):
        return 0.0, 0.0, 0.0, 0.5, -np.inf, -np.inf, -np.inf

    Re = ranking_matrix(E)
    Q = _coranking_matrix(Rx, Re)
    Q = Q[1:, 1:]
    m = len(Q)

    return Q, m


@metric("continuity", paper_reference="zhang2021pydrmetrics", maximize=True)
def continuity(adata: AnnData) -> float:
    Q, m = _fit(adata)
    C = _continuity(Q, m)[_K]
    return float(np.clip(C, 0.0, 1.0))  # in [0, 1]


@metric("co-KNN size", paper_reference="zhang2021pydrmetrics", maximize=True)
def qnn(adata: AnnData) -> float:
    Q, m = _fit(adata)
    QNN = _qnn(Q, m)[_K]
    # normalized in the code to [0, 1]
    return float(np.clip(QNN, 0.0, 1.0))


@metric("co-KNN AUC", paper_reference="zhang2021pydrmetrics", maximize=True)
def qnn_auc(adata: AnnData) -> float:
    Q, m = _fit(adata)
    QNN = _qnn(Q, m)
    AUC = _qnn_auc(QNN)
    return float(np.clip(AUC, 0.5, 1.0))  # in [0.5, 1]


@metric(
    "local continuity meta criterion",
    paper_reference="zhang2021pydrmetrics",
    maximize=True,
)
def lcmc(adata: AnnData) -> float:
    Q, m = _fit(adata)
    QNN = _qnn(Q, m)
    LCMC = _lcmc(QNN, m)[_K]
    return LCMC


@metric("local property", paper_reference="zhang2021pydrmetrics", maximize=True)
def qlocal(adata: AnnData) -> float:
    # according to authors, this is usually preferred to
    # qglobal, because human are more sensitive to nearer neighbors
    Q, m = _fit(adata)
    QNN = _qnn(Q, m)
    LCMC = _lcmc(QNN, m)
    kmax = _kmax(LCMC)
    Qlocal = _q_local(QNN, kmax)
    return Qlocal


@metric("global property", paper_reference="zhang2021pydrmetrics", maximize=True)
def qglobal(adata: AnnData) -> float:
    Q, m = _fit(adata)
    QNN = _qnn(Q, m)
    LCMC = _lcmc(QNN, m)
    kmax = _kmax(LCMC)
    Qglobal = _q_global(QNN, kmax, m)
    return Qglobal
