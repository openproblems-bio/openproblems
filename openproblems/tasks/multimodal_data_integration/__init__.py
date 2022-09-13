from ... import utils
from . import api
from . import datasets
from . import methods
from . import metrics

_task_name = "Multimodal Data Integration"
_task_summary = (
    "Jointly embedding single-cell omics datasets across modalities and batches"
)

DATASETS = utils.get_callable_members(datasets)
METHODS = utils.get_callable_members(methods)
METRICS = utils.get_callable_members(metrics)
