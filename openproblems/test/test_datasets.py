import anndata
import parameterized
import openproblems
from openproblems.test import utils

utils.ignore_numba_warnings()


@parameterized.parameterized.expand(
    [
        (dataset, task, test)
        for task in openproblems.TASKS
        for dataset in task.DATASETS
        for test in [True, False]
    ],
    name_func=utils.name_test,
)
def test_task(dataset, task, test):
    X = dataset(test=test)
    assert isinstance(X, anndata.AnnData)
    assert X.shape[0] > 0
    assert X.shape[1] > 0
    if test:
        assert X.shape[0] < 600
        assert X.shape[1] < 1500
    assert task.checks.check_dataset(X)


@parameterized.parameterized.expand(
    [(dataset,) for task in openproblems.TASKS for dataset in task.DATASETS],
    name_func=utils.name_test,
)
def test_dataset_metadata(dataset):
    assert hasattr(dataset, "metadata")
    for attr in [
        "dataset_name",
    ]:
        assert attr in dataset.metadata
