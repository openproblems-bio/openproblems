from ... import utils
from . import datasets, methods, metrics

_task_name = "Multimodal Data Integration"

DATASETS = utils.get_callable_members(datasets)
METHODS = utils.get_callable_members(methods)
METRICS = utils.get_callable_members(metrics)
