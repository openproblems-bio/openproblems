from .....tools.decorators import dataset
from .....tools.normalize import log_scran_pooling
from ...batch_integration_graph.datasets.immune import immune_batch
from ...batch_integration_graph.datasets.pancreas import pancreas_batch
from typing import Callable


def _convert_dataset_function(
    func: Callable, name: str, image: str = "openproblems-r-base"
) -> Callable:
    metadata = func.metadata.copy()
    metadata["image"] = image

    @dataset(**metadata)
    def converted_func(test=False):
        adata = func(test=test)
        adata = log_scran_pooling(adata)
        adata.X = adata.layers["counts"]
        adata.var_names_make_unique()
        return adata

    converted_func.__name__ = name

    return converted_func


immune_batch_log_scran = _convert_dataset_function(
    immune_batch, "immune_batch_log_scran"
)
pancreas_batch_log_scran = _convert_dataset_function(
    pancreas_batch, "pancreas_batch_log_scran"
)
