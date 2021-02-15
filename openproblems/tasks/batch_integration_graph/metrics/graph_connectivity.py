import numpy as np

from ....tools.decorators import metric
from scIB.metrics import graph_connectivity


@metric(
    metric_name="Graph connectivity",
    maximize=True,
    # image="openproblems-template-image" # only if required
)
def graph_connectivity(adata):
    return graph_connectivity(adata, 'labels')
