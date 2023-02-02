from ....tools.decorators import metric
from ....tools.normalize import log_cp10k


def _distance_correlation(X, X_emb):
    import scipy.spatial
    import scipy.stats

    high_dimensional_distance_vector = scipy.spatial.distance.pdist(X)
    low_dimensional_distance_vector = scipy.spatial.distance.pdist(X_emb)
    return scipy.stats.spearmanr(
        low_dimensional_distance_vector, high_dimensional_distance_vector
    )[0]


@metric(
    metric_name="Distance correlation",
    metric_summary=(
        "Spearman correlation between all pairwise Euclidean "
        "distances in the original and dimension-reduced data"
    ),
    maximize=True,
    paper_reference="schober2018correlation",
)
def distance_correlation(adata, n_svd=200):
    """Calculate the root mean squared error.

    Computes (RMSE) between the full (or processed) data matrix and the
    dimensionally-reduced matrix, invariant to scalar multiplication
    """
    import sklearn.decomposition

    adata = log_cp10k(adata)

    X = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata.X)
    return _distance_correlation(X, adata.obsm["X_emb"])


@metric(
    metric_name="Distance correlation (spectral)",
    metric_summary=(
        "Spearman correlation between all pairwise diffusion "
        "distances in the original and dimension-reduced data"
    ),
    maximize=True,
    paper_reference="coifman2006diffusion",
)
def distance_correlation_spectral(adata, n_comps=200):
    """Calculate the spectral root mean squared error

    Computes (RMSE) between high-dimensional Laplacian eigenmaps on the full (or
    processed) data matrix and the dimensionally-reduced matrix, invariant to scalar
    multiplication
    """
    import numpy as np
    import umap
    import umap.spectral

    adata = log_cp10k(adata)

    n_comps = min(n_comps, min(adata.shape) - 2)

    graph = umap.UMAP(transform_mode="graph").fit_transform(adata.X)
    X = umap.spectral.spectral_layout(
        adata.X, graph, n_comps, random_state=np.random.default_rng()
    )
    return _distance_correlation(X, adata.obsm["X_emb"])
