import parameterized
import unittest
import numpy as np
from scipy import sparse

import openproblems
from openproblems.test import utils, log_level

openproblems.log.setLevel(log_level)
utils.ignore_warnings()


def _dense_data(X):
    if sparse.issparse(X):
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
    class_name_func=utils.name_test,
)
class TestNormalizeX(unittest.TestCase):
    """Test normalization of adata.X."""

    @classmethod
    def setUpClass(cls):
        """Generate and normalize data."""
        cls.adata = utils.data()
        cls.cache_name = cls.normalizer.__name__
        assert utils.assert_finite(cls.adata.X)
        assert cls.cache_name not in cls.adata.layers
        cls.normalizer(cls.adata)

    def test_finite(self):
        """Test that normalized data is finite."""
        assert utils.assert_finite(self.adata.X)

    def test_layers(self):
        """Test that normalized data is cached in adata.layers."""
        assert self.cache_name in self.adata.layers
        utils.assert_array_equal(self.adata.X, self.adata.layers[self.cache_name])

    def test_cache(self):
        """Test that rerunning normalizer loads cached data."""
        # modify normalized data
        self.adata.layers[self.cache_name] = self.adata.layers[self.cache_name].copy()
        cache_data = _dense_data(self.adata.layers[self.cache_name])
        current_data = _dense_data(self.adata.X)

        cache_data += 1
        assert np.all(current_data != _dense_data(self.adata.layers[self.cache_name]))

        # use cached
        self.normalizer(self.adata)
        utils.assert_array_equal(self.adata.X, self.adata.layers[self.cache_name])


@parameterized.parameterized_class(
    ("normalizer"),
    [
        (staticmethod(normalizer),)
        for normalizer in openproblems.utils.get_callable_members(
            openproblems.tools.normalize
        )
    ],
    class_name_func=utils.name_test,
)
class TestNormalizeObsM(unittest.TestCase):
    """Test normalization of adata.X."""

    @classmethod
    def setUpClass(cls):
        """Generate and normalize data."""
        cls.obsm = "test"
        cls.adata = utils.data(obsm=cls.obsm)
        cls.cache_name = "{}_{}".format(cls.obsm, cls.normalizer.__name__)
        assert utils.assert_finite(cls.adata.obsm[cls.obsm])
        assert cls.cache_name not in cls.adata.obsm
        cls.normalizer(cls.adata, obsm=cls.obsm)

    def test_finite(self):
        """Test that normalized data is finite."""
        assert utils.assert_finite(self.adata.obsm[self.obsm])

    def test_layers(self):
        """Test that normalized data is cached in adata.obsm."""
        assert self.cache_name in self.adata.obsm
        utils.assert_array_equal(
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
        self.normalizer(self.adata, obsm=self.obsm)
        utils.assert_array_equal(
            self.adata.obsm[self.obsm], self.adata.obsm[self.cache_name]
        )
