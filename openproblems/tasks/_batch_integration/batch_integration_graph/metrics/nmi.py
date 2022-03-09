from .....tools.decorators import metric


@metric(
    metric_name="NMI",
    maximize=True,
    image="openproblems-python-batch-integration",  # only if required
)
def nmi(adata):
    from scIB.clustering import opt_louvain
    from scIB.metrics import nmi

    opt_louvain(
        adata,
        label_key="labels",
        cluster_key="cluster",
        plot=False,
        inplace=True,
        force=True,
    )
    return nmi(adata, group1="cluster", group2="labels")
