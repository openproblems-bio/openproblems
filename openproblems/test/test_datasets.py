import anndata
import parameterized
import openproblems
from openproblems.test import utils

utils.ignore_numba_warnings()


@parameterized.parameterized(
    [
        (dataset, task, test)
        for task in openproblems.TASKS
        for dataset in task.DATASETS
        for test in [True, False]
    ]
)
def test_task(dataset, task, test):
    X = dataset(test=test)
    assert isinstance(X, anndata.AnnData)
    assert task.checks.check_dataset(X)
    if test:
        assert X.shape[0] < 600
        assert X.shape[1] < 1500
