import numpy as np
import openproblems
import parameterized
import scipy.sparse
import time
import unittest
import utils.asserts
import utils.data
import utils.name


def _dense_data(X):
    if scipy.sparse.issparse(X):
        return X.data
    else:
        return X


@parameterized.parameterized_class(
    ("normalizer"),
    [
        (staticmethod(normalizer),)
        for normalizer in openproblems.utils.get_callable_members(
            openproblems.tools.normalize
        )
    ],
    class_name_func=utils.name.name_test,
)
class TestNormalizeX(unittest.TestCase):
    """Test normalization of adata.X."""

    @classmethod
    def setUpClass(cls):
        """Generate and normalize data."""
        cls.adata = utils.data.data()
        cls.counts = cls.adata.layers["counts"].copy()
        cls.cache_name = cls.normalizer.__name__
        assert utils.asserts.assert_finite(cls.adata.X)
        assert cls.cache_name not in cls.adata.layers
        cls.adata = cls.normalizer(cls.adata)

    def test_not_inplace(self):
        """Test that normalization does not happen inplace."""
        utils.asserts.assert_array_equal(self.adata.layers["counts"], self.counts)

    def test_finite(self):
        """Test that normalized data is finite."""
        assert utils.asserts.assert_finite(self.adata.X)

    def test_layers(self):
        """Test that normalized data is cached in adata.layers."""
        assert self.cache_name in self.adata.layers
        utils.asserts.assert_array_equal(
            self.adata.X, self.adata.layers[self.cache_name]
        )

    def test_cache(self):
        """Test that rerunning normalizer loads cached data."""
        # modify normalized data
        self.adata.layers[self.cache_name] = self.adata.layers[self.cache_name].copy()
        cache_data = _dense_data(self.adata.layers[self.cache_name])
        current_data = _dense_data(self.adata.X)

        cache_data += 1
        assert np.all(current_data != _dense_data(self.adata.layers[self.cache_name]))

        # use cached
        self.adata = self.normalizer(self.adata)
        utils.asserts.assert_array_equal(
            self.adata.X, self.adata.layers[self.cache_name]
        )


@parameterized.parameterized_class(
    ("normalizer"),
    [
        (staticmethod(normalizer),)
        for normalizer in openproblems.utils.get_callable_members(
            openproblems.tools.normalize
        )
    ],
    class_name_func=utils.name.name_test,
)
class TestNormalizeObsM(unittest.TestCase):
    """Test normalization of adata.X."""

    @classmethod
    def setUpClass(cls):
        """Generate and normalize data."""
        cls.obsm = "test"
        cls.adata = utils.data.data(obsm=cls.obsm)
        cls.cache_name = "{}_{}".format(cls.obsm, cls.normalizer.__name__)
        assert utils.asserts.assert_finite(cls.adata.obsm[cls.obsm])
        assert cls.cache_name not in cls.adata.obsm
        cls.adata = cls.normalizer(cls.adata, obsm=cls.obsm)

    def test_finite(self):
        """Test that normalized data is finite."""
        assert utils.asserts.assert_finite(self.adata.obsm[self.obsm])

    def test_layers(self):
        """Test that normalized data is cached in adata.obsm."""
        assert self.cache_name in self.adata.obsm
        utils.asserts.assert_array_equal(
            self.adata.obsm[self.obsm], self.adata.obsm[self.cache_name]
        )

    def test_cache(self):
        """Test that rerunning normalizer loads cached data."""
        # modify normalized data
        self.adata.obsm[self.cache_name] = self.adata.obsm[self.cache_name].copy()
        cache_data = _dense_data(self.adata.obsm[self.cache_name])
        current_data = _dense_data(self.adata.obsm[self.obsm])

        cache_data += 1
        assert np.all(current_data != _dense_data(self.adata.obsm[self.cache_name]))

        # use cached
        self.adata = self.normalizer(self.adata, obsm=self.obsm)
        utils.asserts.assert_array_equal(
            self.adata.obsm[self.obsm], self.adata.obsm[self.cache_name]
        )


def test_profile():
    """Test the profiler."""

    @openproblems.tools.decorators.profile
    def test_fn():
        X = np.random.normal(0, 1, (100, 100))
        time.sleep(2)
        result = X.sum()
        del X
        return result

    result = test_fn()
    for key in ["result", "runtime_s", "memory_mb", "memory_leaked_mb"]:
        assert key in result, key
        assert isinstance(result[key], np.float64 if key == "result" else float), key

    assert result["memory_mb"] >= 0
    assert result["memory_leaked_mb"] >= 0
    assert result["memory_leaked_mb"] < result["memory_mb"] / 100
    assert result["runtime_s"] >= 2
    assert result["runtime_s"] < 3
