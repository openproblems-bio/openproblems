from ....tools.decorators import metric

import anndata
import numpy as np
import scanpy as sc
import sklearn.metrics


@metric(metric_name="Mean-squared error", maximize=False)
def mse(adata):

    test_data = anndata.AnnData(X=adata.obsm["test"], obs=adata.obs, var=adata.var)
    denoised_data = anndata.AnnData(
        X=adata.obsm["denoised"], obs=adata.obs, var=adata.var
    )

    # scaling and transformation
    target_sum = np.median(test_data.X.sum(axis=1))

    sc.pp.normalize_total(test_data, target_sum)
    sc.pp.log1p(test_data)

    sc.pp.normalize_total(denoised_data, target_sum)
    sc.pp.log1p(denoised_data)
    #trying to remove the error in line 26
    mse_res = sklearn.metrics.mean_squared_error(test_data.X.to_array(), denoised_data.X.to_array())
    return mse_res
