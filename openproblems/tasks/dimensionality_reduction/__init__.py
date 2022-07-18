from ... import utils
from . import api
from . import datasets
from . import methods
from . import metrics

_task_name = "Dimensionality reduction for visualisation"
_task_summary = (
    "Reduction of a higher-dimensional feature space to two "
    "dimensions for data visualisation and interpretation"
)

DATASETS = utils.get_callable_members(datasets)
METHODS = utils.get_callable_members(methods)
METRICS = utils.get_callable_members(metrics)
