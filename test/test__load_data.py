import utils

import parameterized
import openproblems

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
    import os

    task = getattr(openproblems.tasks, task_name)
    dataset = getattr(task.datasets, dataset_name)
    adata = dataset(test=True)
    assert not adata.__from_cache__
    adata = dataset(test=True)
    assert adata.__from_cache__
