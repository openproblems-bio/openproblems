from ... import utils
from . import api
from . import datasets
from . import methods
from . import metrics

_task_name = "Label Projection"
_task_summary = "Automated cell type annotation from rich, labeled reference data"

DEFAULT_LAYER = "counts"

DATASETS = utils.get_callable_members(datasets)
METHODS = utils.get_callable_members(methods)
METRICS = utils.get_callable_members(metrics)
