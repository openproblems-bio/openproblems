from .... import utils
from . import api
from . import datasets

_task_name = "Batch integration graph"

DATASETS = utils.get_callable_members(datasets)
