from .....tools.decorators import metric


@metric(
    metric_name="Batch ASW",
    maximize=True,
    image="openproblems-python-batch-integration",  # only if required
)
def silhouette_batch(adata):
    from scib.metrics import silhouette_batch

    sil = silhouette_batch(adata, batch_key="batch", group_key="labels", embed="X_emb")
    return sil
