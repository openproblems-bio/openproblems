import numpy as np

from .....tools.decorators import metric
from scIB.metrics import graph_connectivity


@metric(
    metric_name="Isolated label F1",
    maximize=True,
    # image="openproblems-template-image" # only if required
)
def isolated_labels_f1(adata):
    return isolated_labels(
        adata, label_key="labels", batch_key="batch", cluster=True, verbose=False
    )
