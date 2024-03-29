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
    metric_summary=(
        "Spearman correlation between all pairwise Euclidean distances in the original"
        " and dimension-reduced data"
    ),
    maximize=True,
    paper_reference="schober2018correlation",
)
def distance_correlation(adata, n_svd=500):
    """Calculate the distance correlation

    Computes Spearman correlations between distances on the full (or processed) data
    matrix and the dimensionally-reduced matrix
    """
    import sklearn.decomposition

    adata = log_cp10k(adata)
    X = adata.X
    if n_svd < min(X.shape):
        X = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(X)
    else:
        X = X.toarray()
    return _distance_correlation(X, adata.obsm["X_emb"])


@metric(
    metric_name="Distance correlation (spectral)",
    metric_summary=(
        "Spearman correlation between all pairwise diffusion distances in the original"
        " and dimension-reduced data"
    ),
    maximize=True,
    paper_reference="coifman2006diffusion",
)
def distance_correlation_spectral(adata, n_comps=1000):
    """Calculate the spectral distance correlation

    Computes Spearman correlations between distances on high-dimensional Laplacian
    eigenmaps on the full (or processed) data matrix and the dimensionally-reduced
    matrix
    """
    n_comps = min(n_comps, min(adata.shape) - 2)
    adata_true = diffusion_map(adata.copy(), n_comps=n_comps)
    return _distance_correlation(adata_true.obsm["X_emb"], adata.obsm["X_emb"])
