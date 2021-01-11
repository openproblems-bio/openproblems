import utils

import parameterized
import openproblems

utils.warnings.ignore_warnings()


@parameterized.parameterized.expand(
    [
        (
            task.__name__.split(".")[-1],
            dataset.__name__,
            test,
            dataset.metadata["image"],
            utils.TEMPDIR.name,
        )
        for task in openproblems.TASKS
        for dataset in task.DATASETS
        for test in [True, False]
    ],
    name_func=utils.name.name_test,
)
@utils.docker.docker_test
def test_load_dataset(task_name, dataset_name, test, image, tempdir):
    """Test loading and caching of a dataset."""
    task = getattr(openproblems.tasks, task_name)
    dataset = getattr(task.datasets, dataset_name)
    adata = dataset(test=test)
    adata2 = dataset(test=test)
    assert adata2.shape == adata.shape
    assert adata2.__from_cache__
    utils.cache.save(adata, tempdir, task, dataset, test=test)
