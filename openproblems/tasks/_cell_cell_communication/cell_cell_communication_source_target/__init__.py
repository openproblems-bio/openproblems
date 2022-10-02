from ....utils import get_callable_members
from . import api
from . import datasets
from . import methods
from . import metrics

_task_name = "Cell-Cell Communication Inference (Source-Target)"
_task_summary = "Detect interactions between source and target cell types"

DEFAULT_LAYER = "counts"

DATASETS = get_callable_members(datasets)
METHODS = get_callable_members(methods)
METRICS = get_callable_members(metrics)
