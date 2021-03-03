from .....tools.decorators import metric

import numpy as np


@metric(
    metric_name="Isolated label F1",
    maximize=True,
    image="openproblems-python-batch-integration" # only if required
)
def isolated_labels_f1(adata):
    from scIB.metrics import graph_connectivity
    return isolated_labels(
        adata, label_key="labels", batch_key="batch", cluster=True, verbose=False
    )
