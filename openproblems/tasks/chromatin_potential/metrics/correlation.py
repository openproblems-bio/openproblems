from ....tools.decorators import metric

import numpy as np
import scipy.sparse
import scipy.stats
from sklearn.metrics import mean_squared_error


def _metric(adata, method="pearson"):
    if method == "pearson":
        method = scipy.stats.pearsonr
    elif method == "spearman":
        method = scipy.stats.spearmanr
    elif method == "mse":
        method = mean_squared_error

    metrics = []
    for i in range(adata.shape[0]):
        x = (
            adata.obsm["gene_score"][i].toarray()[0]
            if scipy.sparse.issparse(adata.obsm["gene_score"])
            else adata.obsm["gene_score"][i]
        )
        y = adata.X[i].toarray()[0] if scipy.sparse.issparse(adata.X) else adata.X[i]
        metrics.append(method(x, y)[0])
    metrics = np.array(metrics)
    adata.obs["atac_rna_cor"] = metrics
    return np.median(metrics[~np.isnan(metrics)])


@metric(metric_name="Median Pearson correlation", maximize=True)
def pearson_correlation(adata):
    return _metric(adata)


@metric(metric_name="Median Spearman correlation", maximize=True)
def spearman_correlation(adata):
    return _metric(adata, method="spearman")


@metric(metric_name="Mean squared error", minimize=True)
def mse(adata):
    return _metric(adata, method="mse")
