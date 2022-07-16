from .....tools.decorators import dataset
from .....tools.normalize import log_scran_pooling
from ...batch_integration_graph.datasets.immune import immune_batch
from ...batch_integration_graph.datasets.pancreas import pancreas_batch


@dataset(**(immune_batch.metadata))
def immune_batch_log_scran(test=False):
    adata = immune_batch(test=test)
    return log_scran_pooling(adata)


@dataset(**(pancreas_batch.metadata))
def pancreas_batch_log_scran(test=False):
    adata = pancreas_batch(test=test)
    return log_scran_pooling(adata)
