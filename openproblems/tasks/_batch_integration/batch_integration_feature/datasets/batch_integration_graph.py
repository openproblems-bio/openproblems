from .....tools.decorators import dataset
from .....tools.normalize import log_scran_pooling
from ...batch_integration_graph.datasets.immune import immune_batch
from ...batch_integration_graph.datasets.pancreas import pancreas_batch


def _convert_dataset_function(func):
    metadata = func.metadata.copy()
    metadata["image"] = "openproblems-r-base"

    @dataset(**metadata)
    def converted_func(test=False):
        adata = func(test=test)
        adata = log_scran_pooling(adata)
        adata.X = adata.layers["counts"]
        return adata

    return converted_func


immune_batch_log_scran = _convert_dataset_function(immune_batch)
pancreas_batch_log_scran = _convert_dataset_function(pancreas_batch)
