import anndata
import parameterized
import openproblems


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
