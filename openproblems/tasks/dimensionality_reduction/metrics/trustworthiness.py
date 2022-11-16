from ....tools.decorators import metric
from anndata import AnnData

import numpy as np


@metric(metric_name="trustworthiness", maximize=True)
def trustworthiness(adata: AnnData) -> float:
    from sklearn import manifold

    high_dim, low_dim = adata.X, adata.obsm["X_emb"]

    score = manifold.trustworthiness(
        high_dim, low_dim, n_neighbors=15, metric="euclidean"
    )
    # for large k close to #samples, it's higher than 1.0, e.g 1.0000073552559712
    return float(np.clip(score, 0, 1))
