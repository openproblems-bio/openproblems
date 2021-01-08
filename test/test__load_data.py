import utils

import parameterized
import openproblems
import os

utils.warnings.ignore_warnings()


@parameterized.parameterized.expand(
    [
        (task.__name__.split(".")[-1], dataset.__name__, dataset.metadata["image"])
        for task in openproblems.TASKS
        for dataset in task.DATASETS
    ],
    name_func=utils.name.name_test,
)
@utils.docker.docker_test
def test_load_dataset(task_name, dataset_name, image):
    """Test loading and caching of a dataset."""
    task = getattr(openproblems.tasks, task_name)
    dataset = getattr(task.datasets, dataset_name)
    dataset(test=True)
    assert os.path.isfile(openproblems.data.utils._cache_path(dataset, test=True))
