from ... import utils
from . import api
from . import datasets
from . import methods
from . import metrics

_task_name = "Multimodal Data Integration"
_task_summary = "Alignment of cellular profiles from two different modalities"

DATASETS = utils.get_callable_members(datasets)
METHODS = utils.get_callable_members(methods)
METRICS = utils.get_callable_members(metrics)
