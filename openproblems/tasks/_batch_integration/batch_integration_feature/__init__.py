from .... import utils
from . import api
from . import datasets
from . import methods
from . import metrics

_task_name = "Batch integration embed"
_task_summary = "Feature subtask of batch integration task"

DATASETS = utils.get_callable_members(datasets)
METHODS = utils.get_callable_members(methods)
METRICS = utils.get_callable_members(metrics)
