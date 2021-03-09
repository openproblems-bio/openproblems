from ... import utils
from . import checks
from . import datasets
from . import methods
from . import metrics

# TODO: update
_task_name = "Template Task"

DATASETS = utils.get_callable_members(datasets)
METHODS = utils.get_callable_members(methods)
METRICS = utils.get_callable_members(metrics)
