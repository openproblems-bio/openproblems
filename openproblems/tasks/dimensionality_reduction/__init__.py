from ... import utils
from . import api
from . import datasets
from . import methods
from . import metrics

_task_name = "Dimensionality reduction for visualisation"
_task_summary = (
    "Reduction of high-dimensional datasets to 2D for visualization & interpretation"
)

DEFAULT_LAYER = "log_cp10k"

DATASETS = utils.get_callable_members(datasets)
METHODS = utils.get_callable_members(methods)
METRICS = utils.get_callable_members(metrics)
