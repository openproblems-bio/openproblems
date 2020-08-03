import unittest
import parameterized
import numpy as np

import openproblems
from openproblems.test import utils

utils.ignore_numba_warnings()


@parameterized.parameterized.expand(
    [
        (normalizer,)
        for normalizer in openproblems.utils.get_callable_members(
            openproblems.tools.normalize
        )
    ],
    name_func=utils.name_test,
)
def test_normalize(normalizer):
    adata = utils.data()
    assert normalizer.__name__ not in adata.layers

    # normalize from scratch
    normalizer(adata)
    assert normalizer.__name__ in adata.layers
    np.testing.assert_array_equal(adata.X, adata.layers[normalizer.__name__])

    # modify normalized data
    adata.layers[normalizer.__name__] = adata.layers[normalizer.__name__].copy()
    adata.layers[normalizer.__name__] += 1
    assert np.all(adata.X != adata.layers[normalizer.__name__])

    # use cached
    normalizer(adata)
    np.testing.assert_array_equal(adata.X, adata.layers[normalizer.__name__])
