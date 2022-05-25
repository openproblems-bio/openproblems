from .....tools.decorators import metric


@metric(
    metric_name="Cell Cycle Score",
    maximize=True,
    image="openproblems-python-batch-integration",  # only if required
)
def cc_score(adata):
    from ._utils import _get_split
    from scib.metrics import cell_cycle
    from numpy import nan

    if adata.n_vars<1000:
        return nan

    return cell_cycle(
        *_get_split(adata), "batch", embed="X_emb", organism="human"
    )
