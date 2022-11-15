from .... import utils
from . import api
from . import datasets
from . import methods
from . import metrics

_task_name = "Batch integration feature"
_task_summary = (
    "Removing batch effects while preserving biological variation (feature output)"
)

DEFAULT_LAYER = "log_normalized"

DATASETS = utils.get_callable_members(datasets)
METHODS = utils.get_callable_members(methods)
METRICS = utils.get_callable_members(metrics)
