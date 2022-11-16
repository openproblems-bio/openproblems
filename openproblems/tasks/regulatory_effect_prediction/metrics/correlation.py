from ....tools.decorators import metric

import numpy as np


def _correlation(adata, method="pearson"):
    import scipy.sparse
    import scipy.stats

    if method == "pearson":
        method = scipy.stats.pearsonr
    else:
        method = scipy.stats.spearmanr

    cors = []
    for i in range(adata.shape[0]):
        x = (
            adata.obsm["gene_score"][i].toarray()[0]
            if scipy.sparse.issparse(adata.obsm["gene_score"])
            else adata.obsm["gene_score"][i]
        )
        y = adata.X[i].toarray()[0] if scipy.sparse.issparse(adata.X) else adata.X[i]
        cors.append(method(x, y)[0])
    cors = np.array(cors)
    adata.obs["atac_rna_cor"] = cors
    return np.median(cors[~np.isnan(cors)])


@metric(metric_name="Median Pearson correlation", maximize=True)
def pearson_correlation(adata):
    return _correlation(adata)


@metric(metric_name="Median Spearman correlation", maximize=True)
def spearman_correlation(adata):
    return _correlation(adata, method="spearman")
