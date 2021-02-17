import numpy as np

from .....tools.decorators import metric
import scIB.metrics


@metric(
    metric_name="Graph connectivity",
    maximize=True,
    # image="openproblems-template-image" # only if required
)
def graph_connectivity(adata):
    return scIB.metrics.graph_connectivity(adata, "labels")
