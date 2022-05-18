from .....tools.decorators import metric


@metric(
    metric_name="Isolated label Silhouette",
    maximize=True,
    image="openproblems-python-batch-integration",  # only if required
)
def isolated_labels_sil(adata):
    from scib.metrics import isolated_labels

    return isolated_labels(
        adata,
        label_key="labels",
        batch_key="batch",
        embed="X_emb",
        cluster=False,
        verbose=False,
    )
