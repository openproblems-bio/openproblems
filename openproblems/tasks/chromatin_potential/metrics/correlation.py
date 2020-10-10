
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.sparse import issparse
import numpy as np
from ....tools.decorators import metric


@metric(metric_name="Correlation between RNA and ATAC", maximize=True)
def correlation(adata, method='pearson'):
    if method == 'pearson':
        method = pearsonr
    else:
        method = spearmanr

    cors = []
    for i in range(adata.shape[0]):
        x = adata.obsm['gene_score'][i].toarray()[0] if issparse(adata.obsm['gene_score']) else adata.obsm['gene_score'][i]
        y = adata.X[i].toarray()[0] if issparse(adata.X) else adata.X[i]
        cors.append(method(x, y)[0])
    cors = np.array(cors)
    return cors[~np.isnan(cors)].mean()
