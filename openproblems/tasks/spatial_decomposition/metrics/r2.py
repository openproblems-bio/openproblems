from ....tools.decorators import metric


@metric(
    metric_name="r2",
    metric_summary=(
        'R2, also known as the "coefficient of determination", reports the fraction of'
        " the true proportion values' variance that can be explained by the predicted"
        " proportion values. The **best score**, and upper bound, is 1.0. There is no"
        " fixed lower bound for the metric. The _uniform/non-weighted average_ across"
        " all cell types/states is used to summarize performance."
    ),
    maximize=True,
    paper_reference="miles2005rsquared",
)
def r2(adata):
    import sklearn.metrics

    prop_true = adata.obsm["proportions_true"]
    prop_pred = adata.obsm["proportions_pred"]

    r2_score = sklearn.metrics.r2_score(
        prop_true, prop_pred, sample_weight=None, multioutput="uniform_average"
    )
    return r2_score
