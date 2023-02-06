from anndata import AnnData
from openproblems.tools.decorators import metric
from typing import Optional

import numpy as np

__license__ = """
BSD 3-Clause License

Copyright (c) 2017, Leland McInnes
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
_K = 30  # number of neighbors
_SEED = 42


def _calculate_radii(
    X: np.ndarray, n_neighbors: int = 30, random_state: Optional[int] = None
) -> np.ndarray:
    from umap.umap_ import fuzzy_simplicial_set
    from umap.umap_ import nearest_neighbors

    # directly taken from: https://github.com/lmcinnes/umap/blob/
    # 317ce81dc64aec9e279aa1374ac809d9ced236f6/umap/umap_.py#L1190-L1243
    knn_indices, knn_dists, _ = nearest_neighbors(
        X,
        n_neighbors,
        "euclidean",
        {},
        False,
        random_state,
        verbose=False,
    )

    emb_graph, _, _, emb_dists = fuzzy_simplicial_set(
        X,
        n_neighbors,
        random_state,
        "euclidean",
        {},
        knn_indices,
        knn_dists,
        verbose=False,
        return_dists=True,
    )

    emb_graph = emb_graph.tocoo()
    emb_graph.sum_duplicates()
    emb_graph.eliminate_zeros()

    n_vertices = emb_graph.shape[1]

    mu_sum = np.zeros(n_vertices, dtype=np.float32)
    re = np.zeros(n_vertices, dtype=np.float32)

    head = emb_graph.row
    tail = emb_graph.col
    for i in range(len(head)):
        j = head[i]
        k = tail[i]

        D = emb_dists[j, k]
        mu = emb_graph.data[i]

        re[j] += mu * D
        re[k] += mu * D
        mu_sum[j] += mu
        mu_sum[k] += mu

    epsilon = 1e-8
    return np.log(epsilon + (re / mu_sum))


@metric(
    "density preservation",
    paper_reference="narayan2021assessing",
    maximize=True,
)
def density_preservation(adata: AnnData) -> float:
    from scipy.sparse import issparse
    from scipy.stats import pearsonr

    emb = adata.obsm["X_emb"]

    high_dim = adata.X.A if issparse(adata.X) else adata.X
    ro = _calculate_radii(high_dim, n_neighbors=_K, random_state=_SEED)
    # in principle, we could just call _calculate_radii(high_dim, ...)
    # this is made sure that the test pass (otherwise, there was .02 difference in corr)
    re = _calculate_radii(emb, n_neighbors=_K, random_state=_SEED)

    return pearsonr(ro, re)[0]
