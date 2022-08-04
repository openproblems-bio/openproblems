from ...utils import get_callable_members
from . import api
from . import datasets
from . import methods
from . import metrics

_task_name = "Cell-Cell Communication Inference"
_task_summary = "Detect intercellular communication events between cell types"

DATASETS = get_callable_members(datasets)
METHODS = get_callable_members(methods)
METRICS = get_callable_members(metrics)
