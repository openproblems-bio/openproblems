from ... import utils
from . import api
from . import datasets
from . import methods
from . import metrics

_task_name = "Denoising"
_task_summary = "Denoising UMI counts by Molecular Cross Validation"

DATASETS = utils.get_callable_members(datasets)
METHODS = utils.get_callable_members(methods)
METRICS = utils.get_callable_members(metrics)
