from ....tools.decorators import metric

import sklearn.metrics
import anndata
import numpy as np


@metric(metric_name="R2", maximize=True)
def r2(adata):

    prop_true = adata.obsm["proportions_true"].values
    prop_pred = adata.obsm["proportions_pred"].values

    r2_score = sklearn.metrics.r2_score(prop_true,
                                        prop_pred,
                                        sample_weight=None,
                                        multioutput='uniform_average')
    return r2_score
