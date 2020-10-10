import numpy as np
from sklearn.metrics import mean_squared_error
import scipy as sc
from scipy import AnnData

from ....tools.decorators import metric

@metric(metric_name="Mean-squared error", maximize=False)
def mse(adata):

    test_data = AnnData(X=adata.obsm["test"], obs=adata.obs, var=adata.var)
    denoised_data = AnnData(X=adata.obsm["denoised"], obs=adata.obs, var=adata.var)

    # scaling and transformation
    target_sum = np.median(test_data.X.sum(axis=1))

    sc.pp.normalize_total(test_data, target_sum)
    sc.pp.log1p(test_data)

    sc.pp.normalize_total(denoised_data, target_sum)
    sc.pp.log1p(denoised_data)

    error = mean_squared_error(test_data.X, denoised_data.X)
    return error
