from ....tools.decorators import metric


@metric(metric_name="r2", maximize=True, paper_reference="miles2005rsquared")
def r2(adata):
    import sklearn.metrics

    prop_true = adata.obsm["proportions_true"]
    prop_pred = adata.obsm["proportions_pred"]

    r2_score = sklearn.metrics.r2_score(
        prop_true, prop_pred, sample_weight=None, multioutput="uniform_average"
    )
    return r2_score
