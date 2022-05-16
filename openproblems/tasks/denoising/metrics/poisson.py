from ....tools.decorators import metric

import scprep


@metric(metric_name="Poisson loss", maximize=False, image="openproblems-python-extras")
def poisson(adata):
    from molecular_cross_validation.mcv_sweep import poisson_nll_loss

    test_data = adata.obsm["test"]
    denoised_data = adata.obsm["denoised"]

    # scaling
    initial_sum = adata.obsm["train"].sum()
    target_sum = test_data.sum()
    denoised_data = denoised_data * target_sum / initial_sum

    error = poisson_nll_loss(scprep.utils.toarray(test_data), denoised_data)
    return error
