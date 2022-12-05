from ....tools.decorators import metric


def _rmse(X, X_emb):
    import scipy.optimize
    import scipy.spatial

    high_dimensional_distance_vector = scipy.spatial.distance.pdist(X)
    low_dimensional_distance_vector = scipy.spatial.distance.pdist(X_emb)
    scale, rmse = scipy.optimize.nnls(
        low_dimensional_distance_vector[:, None], high_dimensional_distance_vector
    )
    return rmse


@metric(metric_name="RMSE", maximize=False)
def rmse(adata, n_svd=200):
    """Calculate the root mean squared error.

    Computes (RMSE) between the full (or processed) data matrix and the
    dimensionally-reduced matrix, invariant to scalar multiplication
    """
    import sklearn.decomposition

    X = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.X)
    return _rmse(X, adata.obsm["X_emb"])
