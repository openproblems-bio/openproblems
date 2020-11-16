import pandas as pd
import anndata
import parameterized
import openproblems
import warnings

import packaging.version
import scipy.sparse
from openproblems.test import utils

utils.ignore_warnings()


def _assert_not_bytes(X):
    if isinstance(X, pd.Series):
        assert not X.apply(lambda x: isinstance(x, bytes)).any()
    elif isinstance(X, pd.Index):
        return _assert_not_bytes(X.to_series())
    elif isinstance(X, pd.DataFrame):
        assert _assert_not_bytes(X.index)
        assert X.apply(_assert_not_bytes).all()
    elif isinstance(X, dict):
        for v in X.values():
            assert _assert_not_bytes(v)
    else:
        pass
    return True


@parameterized.parameterized.expand(
    [
        (dataset, task, test)
        for task in openproblems.TASKS
        for dataset in task.DATASETS
        for test in [True, False]
    ],
    name_func=utils.name_test,
)
def test_dataset(dataset, task, test):
    """Test dataset loading."""
    adata = dataset(test=test)
    assert isinstance(adata, anndata.AnnData)
    assert adata.shape[0] > 0
    assert adata.shape[1] > 0

    is_sparse = scipy.sparse.issparse(adata.X)
    assert is_sparse or packaging.version.parse(
        openproblems.__version__
    ) < packaging.version.parse("1.0")
    if not is_sparse:
        warnings.warn(
            "{}-{}: Dense data will raise an error in openproblems v1.0".format(
                task.__name__.split(".")[-1], dataset.__name__
            ),
            FutureWarning,
        )

    assert _assert_not_bytes(adata.obs)
    assert _assert_not_bytes(adata.var)
    assert _assert_not_bytes(adata.uns)

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
    """Test for existence of dataset metadata."""
    assert hasattr(dataset, "metadata")
    for attr in [
        "dataset_name",
    ]:
        assert attr in dataset.metadata
