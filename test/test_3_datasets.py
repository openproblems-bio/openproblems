import utils.warnings  # noqa: F401

# isort: split
import anndata
import openproblems
import openproblems.utils
import pandas as pd
import parameterized
import pytest
import scipy.sparse
import unittest
import utils
import utils.asserts
import utils.cache
import utils.git
import utils.name

DATASET_SUMMARY_MAXLEN = 80

pytestmark = pytest.mark.skipif(
    len(utils.git.list_modified_tasks()) == 0, reason="No tasks have been modified"
)


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


@parameterized.parameterized_class(
    ("dataset", "task", "test", "tempdir"),
    [
        (staticmethod(dataset), task, test, utils.TEMPDIR.name)
        for task in utils.git.list_modified_tasks()
        for dataset in task.DATASETS
        for test in [True]
    ],
    class_name_func=utils.name.name_test,
)
class TestDataset(unittest.TestCase):
    """Test cached dataset characteristics."""

    @classmethod
    def setUpClass(cls):
        """Load data."""
        try:
            cls.adata = utils.cache.load(
                cls.tempdir,
                cls.task,
                cls.dataset,
                test=cls.test,
                dependency="test_load_dataset",
            )
        except AssertionError as e:
            if str(e) == "Intermediate file missing. Did test_load_dataset fail?":
                pytest.skip("Dataset not loaded successfully")
            else:
                raise

    @classmethod
    def tearDownClass(cls):
        """Remove data."""
        utils.cache.delete(
            cls.tempdir,
            cls.task,
            cls.dataset,
            test=cls.test,
        )

    def test_adata_class(self):
        """Ensure output is AnnData."""
        assert isinstance(self.adata, anndata.AnnData)

    def test_adata_shape(self):
        """Ensure output is of appropriate size."""
        assert self.adata.shape[0] > 0
        assert self.adata.shape[1] > 0
        if self.test:
            assert self.adata.shape[0] <= 600
            assert self.adata.shape[1] <= 1500

    def test_sparse(self):
        """Ensure output is sparse."""
        if not scipy.sparse.issparse(self.adata.X):
            openproblems.utils.future_warning(
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
        assert self.task.api.check_dataset(self.adata)

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
            adata = normalizer(adata)
            utils.asserts.assert_finite(adata.X)

    def test_metadata(self):
        """Test for existence of dataset metadata."""
        assert hasattr(self.dataset, "metadata")
        for attr in ["dataset_name", "data_url", "dataset_summary", "image"]:
            assert attr in self.dataset.metadata
            assert self.dataset.metadata[attr] is not None

        assert isinstance(self.dataset.metadata["dataset_name"], str)
        assert isinstance(self.dataset.metadata["image"], str)
        assert self.dataset.metadata["image"].startswith("openproblems")
        assert isinstance(self.dataset.metadata["dataset_summary"], str)
        assert len(self.dataset.metadata["dataset_summary"]) < DATASET_SUMMARY_MAXLEN
        assert isinstance(self.dataset.metadata["data_url"], str)
        assert utils.asserts.assert_url_accessible(self.dataset.metadata["data_url"])
