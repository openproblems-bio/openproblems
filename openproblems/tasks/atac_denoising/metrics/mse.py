from ....tools.decorators import metric
from sklearn.metrics import mean_squared_error


@metric(metric_name="Mean-squared error", maximize=False)
def mse(adata):
    return mean_squared_error(
        adata.obsm["mode2_denoised"].todense(),
        adata.obsm["mode2"].todense(),
    )
