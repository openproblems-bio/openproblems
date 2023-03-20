from ....tools.decorators import metric


@metric(
    metric_name="Poisson loss",
    metric_summary=(
        "The Poisson log likelihood of observing the true counts of the test dataset"
        " given the distribution given in the denoised dataset."
    ),
    paper_reference="batson2019molecular",
    maximize=False,
    image="openproblems-python-pytorch",
)
def poisson(adata):
    from molecular_cross_validation.mcv_sweep import poisson_nll_loss

    import scprep

    test_data = adata.obsm["test"]
    denoised_data = adata.obsm["denoised"]

    # scaling
    initial_sum = adata.obsm["train"].sum()
    target_sum = test_data.sum()
    denoised_data = denoised_data * target_sum / initial_sum

    error = poisson_nll_loss(scprep.utils.toarray(test_data), denoised_data)
    return error
