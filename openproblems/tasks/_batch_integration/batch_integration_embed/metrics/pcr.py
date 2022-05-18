from .....tools.decorators import metric


@metric(
    metric_name="PC Regression",
    maximize=True,
    image="openproblems-python-batch-integration",  # only if required
)
def pcr(adata):
    from ._utils import _get_split
    from scib.metrics import pcr_comparison

    return pcr_comparison(*_get_split(adata), "batch", embed="X_emb")
