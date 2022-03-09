from ....tools.decorators import metric
from sklearn.metrics import mean_squared_error

import numpy as np
import scipy.sparse
import scipy.stats


def _metric(adata, method="pearson"):
    if method == "pearson":
        method = scipy.stats.pearsonr
    elif method == "spearman":
        method = scipy.stats.spearmanr
    elif method == "mse":
        method = mean_squared_error

    results = []
    for i in range(adata.shape[0]):
        x = (
            adata.obsm["gene_score"][i].toarray()[0]
            if scipy.sparse.issparse(adata.obsm["gene_score"])
            else adata.obsm["gene_score"][i]
        )
        y = adata.X[i].toarray()[0] if scipy.sparse.issparse(adata.X) else adata.X[i]
        result = method(x, y)
        if not isinstance(result, np.float64):
            results.append(result[0])
        else:
            results.append(result)
    results = np.array(results)
    return np.median(results[~np.isnan(results)])


@metric(metric_name="Median Pearson correlation", maximize=True)
def pearson_correlation(adata):
    return _metric(adata, method="pearson")


@metric(metric_name="Median Spearman correlation", maximize=True)
def spearman_correlation(adata):
    return _metric(adata, method="spearman")


@metric(metric_name="Mean squared error", maximize=False)
def mse(adata):
    return _metric(adata, method="mse")
