from ....tools.decorators import metric
from ....tools.normalize import log_cp10k
from ..methods.diffusion_map import diffusion_map


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
    maximize=True,
    paper_reference="schober2018correlation",
)
def distance_correlation(adata, n_svd=1000):
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
    maximize=True,
    paper_reference="coifman2006diffusion",
)
def distance_correlation_spectral(adata, n_comps=1000):
    """Calculate the spectral root mean squared error

    Computes (RMSE) between high-dimensional Laplacian eigenmaps on the full (or
    processed) data matrix and the dimensionally-reduced matrix, invariant to scalar
    multiplication
    """
    n_comps = min(n_comps, min(adata.shape) - 2)
    adata_true = diffusion_map(adata.copy(), n_comps=n_comps)
    return _distance_correlation(adata_true.obsm["X_emb"], adata.obsm["X_emb"])
