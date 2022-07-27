from ... import utils
from . import api
from . import datasets
from . import methods
from . import metrics

_task_name = "Cell-Cell Communication Inference"
_task_summary = (
    "Evaluate the methods' ability to detect assumed truth intercellular events"
)

DATASETS = utils.get_callable_members(datasets)
METHODS = utils.get_callable_members(methods)
METRICS = utils.get_callable_members(metrics)
