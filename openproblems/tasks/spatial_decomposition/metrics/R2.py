from ....tools.decorators import metric

import sklearn.metrics


@metric(metric_name="r2", maximize=True)
def r2(adata):

    prop_true = adata.obsm["proportions_true"]
    prop_pred = adata.obsm["proportions_pred"]

    r2_score = sklearn.metrics.r2_score(
        prop_true, prop_pred, sample_weight=None, multioutput="uniform_average"
    )
    return r2_score
