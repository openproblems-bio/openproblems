import anndata
import parameterized
import openproblems


@parameterized.parameterized(
    [(dataset, task) for task in openproblems.TASKS for dataset in task.DATASETS]
)
def test_task(dataset, task):
    X = dataset()
    assert isinstance(X, anndata.AnnData)
    assert task.checks.check_dataset(X)
