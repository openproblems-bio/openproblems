from .....tools.decorators import metric


@metric(
    metric_name="Isolated label F1",
    maximize=True,
    image="openproblems-python-batch-integration",  # only if required
)
def isolated_labels_f1(adata):
    from scIB.metrics import isolated_labels

    return isolated_labels(
        adata,
        label_key="labels",
        batch_key="batch",
        embed="X_pca",
        cluster=True,
        verbose=False,
    )
