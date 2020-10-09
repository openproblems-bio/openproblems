
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from ....tools.decorators import metric


@metric(metric_name="Correlation between RNA and ATAC", maximize=True)
def correlation(x, y, method='pearson'):
    if method == 'pearson':
        return pearsonr(x, y)[0]
    else:
        return spearmanr(x, y)[0]
