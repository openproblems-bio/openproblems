import utils

import pandas as pd
import anndata
import parameterized
import openproblems
import unittest

import scipy.sparse

utils.warnings.ignore_warnings()


def _assert_not_bytes(X):
    if isinstance(X, pd.Series):
        if pd.api.types.is_categorical_dtype(X):
            return _assert_not_bytes(X.dtype.categories)
        else:
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


@parameterized.parameterized_class(
    ("dataset", "task", "test"),
    [
        (staticmethod(dataset), task, test)
        for task in openproblems.TASKS
        for dataset in task.DATASETS
        for test in [True, False]
    ],
    class_name_func=utils.name.name_test,
)
class TestDataset(unittest.TestCase):
    """Test cached dataset characteristics."""

    @classmethod
    def setUpClass(cls):
        """Load data."""
        cls.adata = cls.dataset(test=cls.test)

    def test_adata_class(self):
        """Ensure output is AnnData."""
        assert isinstance(self.adata, anndata.AnnData)

    def test_adata_shape(self):
        """Ensure output is of appropriate size."""
        assert self.adata.shape[0] > 0
        assert self.adata.shape[1] > 0
        if self.test:
            assert self.adata.shape[0] < 600
            assert self.adata.shape[1] < 1500

    def test_sparse(self):
        """Ensure output is sparse."""
        if not scipy.sparse.issparse(self.adata.X):
            utils.warnings.future_warning(
                "{}-{}: self.adata.X is loaded as dense.".format(
                    self.task.__name__.split(".")[-1], self.dataset.__name__
                ),
                error_version="1.0",
                error_category=TypeError,
            )

    def test_not_bytes(self):
        """Ensure output does not contain byte strings."""
        assert _assert_not_bytes(self.adata.obs)
        assert _assert_not_bytes(self.adata.var)
        assert _assert_not_bytes(self.adata.uns)

    def test_not_empty(self):
        """Ensure output does not have empty rows or columns."""
        assert self.adata.X.sum(axis=0).min() > 0
        assert self.adata.X.sum(axis=1).min() > 0

    def test_task_checks(self):
        """Run task-specific tests."""
        assert self.task.checks.check_dataset(self.adata)

    def test_cache(self):
        """Test that AnnData files are written to disk."""
        adata = self.dataset(test=self.test)
        assert adata.__from_cache__
        assert adata.shape == self.adata.shape

    @parameterized.parameterized.expand(
        [
            (normalizer,)
            for normalizer in openproblems.utils.get_callable_members(
                openproblems.tools.normalize
            )
        ]
    )
    def test_normalize(self, normalizer):
        """Test that normalizations can be safely applied."""
        if self.test:
            adata = self.adata.copy()
            normalizer(adata)
            utils.asserts.assert_finite(adata.X)

    def test_metadata(self):
        """Test for existence of dataset metadata."""
        assert hasattr(self.dataset, "metadata")
        for attr in [
            "dataset_name",
        ]:
            assert attr in self.dataset.metadata
