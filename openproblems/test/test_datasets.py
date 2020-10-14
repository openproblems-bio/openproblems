import anndata
import parameterized
import openproblems
import warnings
from scipy import sparse
from openproblems.test import utils

utils.ignore_warnings()


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
    adata = dataset(test=test)
    assert isinstance(adata, anndata.AnnData)
    assert adata.shape[0] > 0
    assert adata.shape[1] > 0
    if not sparse.issparse(adata.X):
        warnings.warn(
            "{}-{}: Dense data will raise an error in openproblems v1.0".format(
                task.__name__.split(".")[-1], dataset.__name__
            ),
            FutureWarning,
        )
    # assert sparse.issparse(adata.X)
    assert adata.X.sum(axis=0).min() > 0
    assert adata.X.sum(axis=1).min() > 0
    if test:
        assert adata.shape[0] < 600
        assert adata.shape[1] < 1500
    assert task.checks.check_dataset(adata)


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
