from ....tools.decorators import metric
from ....tools.normalize import log_cpm


def _rmse(X, X_emb):
    import scipy.optimize
    import scipy.spatial

    high_dimensional_distance_vector = scipy.spatial.distance.pdist(X)
    low_dimensional_distance_vector = scipy.spatial.distance.pdist(X_emb)
    _, rmse = scipy.optimize.nnls(
        low_dimensional_distance_vector[:, None], high_dimensional_distance_vector
    )
    return rmse


@metric(
    metric_name="RMSE",
    maximize=False,
    paper_reference="kruskal1964mds",
)
def rmse(adata, n_svd=200):
    """Calculate the root mean squared error.

    Computes (RMSE) between the full (or processed) data matrix and the
    dimensionally-reduced matrix, invariant to scalar multiplication
    """
    import sklearn.decomposition

    adata = log_cpm(adata)

    X = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.X)
    return _rmse(X, adata.obsm["X_emb"])


@metric(
    metric_name="RMSE (spectral)",
    maximize=False,
    paper_reference="coifman2006diffusion",
)
def rmse_spectral(adata, n_comps=200):
    """Calculate the spectral root mean squared error

    Computes (RMSE) between high-dimensional Laplacian eigenmaps on the full (or
    processed) data matrix and the dimensionally-reduced matrix, invariant to scalar
    multiplication
    """
    import numpy as np
    import umap
    import umap.spectral

    adata = log_cpm(adata)

    n_comps = min(n_comps, min(adata.shape) - 2)

    graph = umap.UMAP(transform_mode="graph").fit_transform(adata.X)
    X = umap.spectral.spectral_layout(
        adata.X, graph, n_comps, random_state=np.random.default_rng()
    )
    return _rmse(X, adata.obsm["X_emb"])
