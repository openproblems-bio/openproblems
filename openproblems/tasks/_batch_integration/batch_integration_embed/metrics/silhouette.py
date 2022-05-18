from .....tools.decorators import metric


@metric(
    metric_name="Silhouette",
    maximize=True,
    image="openproblems-python-batch-integration",  # only if required
)
def silhouette(adata):
    from scib.metrics import silhouette

    return silhouette(adata, group_key="labels", embed="X_emb")
