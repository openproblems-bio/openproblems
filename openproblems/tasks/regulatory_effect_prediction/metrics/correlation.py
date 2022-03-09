from ....tools.decorators import metric
from sklearn.metrics import mean_squared_error

import numpy as np
import scipy.sparse
import scipy.stats


def _metric(adata, method="pearson", modality="global"):
    if modality == "median":
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
            y = (
                adata.X[i].toarray()[0]
                if scipy.sparse.issparse(adata.X)
                else adata.X[i]
            )
            result = method(x, y)
            if not isinstance(result, np.float64):
                results.append(result[0])
            else:
                results.append(result)
        results = np.array(results)
        return np.median(results[~np.isnan(results)])
    elif modality == "global":
        if method == "pearson":
            return _global_pearson(adata.obsm["gene_score"], adata.X)
        elif method == "spearman":
            return _global_spearman(adata.obsm["gene_score"], adata.X)
        elif method == "mse":
            return mean_squared_error(adata.obsm["gene_score"], adata.X)


def _global_pearson(x, y):
    numerator = np.mean((x - x.mean()) * (y - y.mean()))
    denominator = x.std() * y.std()
    if denominator == 0:
        return 0
    else:
        result = numerator / denominator
        return result


def _global_spearman(x, y):
    x = scipy.stats.rankdata(x)
    y = scipy.stats.rankdata(y)
    numerator = np.mean((x - x.mean()) * (y - y.mean()))
    denominator = x.std() * y.std()
    if denominator == 0:
        return 0
    else:
        result = numerator / denominator
        return result


@metric(metric_name="Global Pearson correlation", maximize=True)
def global_pearson_correlation(adata):
    return _metric(adata, method="pearson", modality="global")


@metric(metric_name="Global Spearman correlation", maximize=True)
def global_spearman_correlation(adata):
    return _metric(adata, method="spearman", modality="global")


@metric(metric_name="Global Mean squared error", maximize=True)
def global_mse(adata):
    return _metric(adata, method="mse", modality="global")


@metric(metric_name="Median Pearson correlation", maximize=True)
def median_pearson_correlation(adata):
    return _metric(adata, method="pearson", modality="median")


@metric(metric_name="Median Spearman correlation", maximize=True)
def median_spearman_correlation(adata):
    return _metric(adata, method="spearman", modality="median")


@metric(metric_name="Mean squared error", maximize=False)
def median_mse(adata):
    return _metric(adata, method="mse", modality="median")
