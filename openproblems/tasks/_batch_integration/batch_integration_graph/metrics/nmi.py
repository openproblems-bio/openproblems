from .....tools.decorators import metric

import numpy as np


@metric(
    metric_name="NMI",
    maximize=True,
    image="openproblems-python-batch-integration",  # only if required
)
def nmi(adata):
    from scIB.clustering import opt_louvain
    from scIB.metrics import nmi

    res_max, nmi_max, nmi_all = opt_louvain(
        adata,
        label_key="labels",
        cluster_key="cluster",
        function=nmi,
        plot=False,
        verbose=verbose,
        inplace=True,
        force=True,
    )
    return nmi(adata, group1="cluster", group2="labels")
