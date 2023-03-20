from ....tools.decorators import metric


@metric(
    metric_name="Mean-squared error",
    metric_summary=(
        "The mean squared error between the denoised counts of the training dataset and"
        " the true counts of the test dataset after reweighting by the train/test"
        " ratio."
    ),
    paper_reference="batson2019molecular",
    maximize=False,
)
def mse(adata):
    import anndata
    import scanpy as sc
    import scprep
    import sklearn.metrics

    test_data = anndata.AnnData(X=adata.obsm["test"], obs=adata.obs, var=adata.var)
    denoised_data = anndata.AnnData(
        X=adata.obsm["denoised"], obs=adata.obs, var=adata.var
    )

    # scaling and transformation
    target_sum = 10000

    sc.pp.normalize_total(test_data, target_sum)
    sc.pp.log1p(test_data)

    sc.pp.normalize_total(denoised_data, target_sum)
    sc.pp.log1p(denoised_data)

    error = sklearn.metrics.mean_squared_error(
        scprep.utils.toarray(test_data.X), denoised_data.X
    )
    return error
