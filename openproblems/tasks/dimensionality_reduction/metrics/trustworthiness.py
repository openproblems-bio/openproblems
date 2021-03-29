from ....tools.decorators import metric
from anndata import AnnData
from scipy.sparse import issparse
from sklearn import manifold

import numpy as np


@metric(metric_name="trustworthiness", maximize=True)
def trustworthiness(adata: AnnData) -> float:
    high_dim, low_dim = adata.X, adata.obsm["X_emb"]

    data = high_dim.data if issparse(high_dim) else high_dim
    if np.any(~np.isfinite(data)):
        # fault of the data
        return float(np.nan)

    if np.any(~np.isfinite(low_dim)):
        # method's fault
        return 0.0

    score = manifold.trustworthiness(
        high_dim, low_dim, n_neighbors=15, metric="euclidean"
    )
    # for large k close to #samples, it's higher than 1.0, e.g 1.0000073552559712
    score = float(np.clip(score, 0, 1))
    adata.uns["trustworthiness_score"] = score

    return score
