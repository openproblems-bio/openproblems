from ... import utils
from . import datasets, methods, metrics, checks

_task_name = "Label Projection"

DATASETS = utils.get_members(datasets)
METHODS = utils.get_members(methods)
METRICS = utils.get_members(metrics)
