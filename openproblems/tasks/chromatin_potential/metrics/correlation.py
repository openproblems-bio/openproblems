
from scipy.stats import pearsonr
from scipy.stats import spearmanr
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
        cors.append(method(adata.obsm['gene_score'][i].toarray()[0], adata.X[i].toarray()[0]))
    return np.median(cors)
