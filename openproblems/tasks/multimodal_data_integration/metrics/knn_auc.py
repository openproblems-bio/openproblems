import numpy as np
from sklearn.neighbors import NearestNeighbors


def knn_auc(adata, n_neighbors=100, n_svd=100):
    if min(adata.X.shape) <= n_svd:
        n_svd = min(adata.X.shape) - 1
    X_pca = TruncatedSVD(n_svd).fit_transform(adata.X)
    _, indices_true = (
        NearestNeighbors(n_neighbors=n_neighbors).fit(X_pca).kneighbors(X_pca)
    )
    _, indices_pred = (
        NearestNeighbors(n_neighbors=n_neighbors)
        .fit(adata.obsm["aligned"])
        .kneighbors(adata.obsm["mode2_aligned"])
    )
    neighbors_match = np.zeros(n_neighbors)
    for i in range(adata.shape[0]):
        neighbors_match += [
            np.sum(np.isin(indices_pred[i, : j + 1], indices_true[i, : j + 1]))
            for j in range(n_neighbors)
        ]

    neighbors_match_curve = neighbors_match / (
        np.arange(1, n_neighbors + 1) * adata.shape[0]
    )
    area_under_curve = np.mean(neighbors_match_curve)
    return area_under_curve
