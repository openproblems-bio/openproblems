import parameterized
import numpy as np
from scipy import sparse

import openproblems
from openproblems.test import utils

utils.ignore_warnings()


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
    utils.assert_array_equal(adata.X, adata.layers[normalizer.__name__])

    # modify normalized data
    adata.layers[normalizer.__name__] = adata.layers[normalizer.__name__].copy()
    if sparse.issparse(adata.layers[normalizer.__name__]):
        adata.layers[normalizer.__name__].data += 1
        assert np.all(adata.X.data != adata.layers[normalizer.__name__].data)
    else:
        adata.layers[normalizer.__name__] += 1
        assert np.all(adata.X != adata.layers[normalizer.__name__])

    # use cached
    normalizer(adata)
    utils.assert_array_equal(adata.X, adata.layers[normalizer.__name__])


@parameterized.parameterized.expand(
    [
        (normalizer,)
        for normalizer in openproblems.utils.get_callable_members(
            openproblems.tools.normalize
        )
    ],
    name_func=utils.name_test,
)
def test_normalize_obsm(normalizer, obsm="test"):
    adata = utils.data(obsm=obsm)
    cache_name = "{}_{}".format(obsm, normalizer.__name__)
    assert cache_name not in adata.obsm

    # normalize from scratch
    normalizer(adata, obsm=obsm)
    assert cache_name in adata.obsm
    utils.assert_array_equal(adata.obsm[obsm], adata.obsm[cache_name])

    # modify normalized data
    adata.obsm[cache_name] = adata.obsm[cache_name].copy()
    if sparse.issparse(adata.obsm[cache_name]):
        adata.obsm[cache_name].data += 1
        assert np.all(adata.obsm[obsm].data != adata.obsm[cache_name].data)
    else:
        adata.obsm[cache_name] += 1
        assert np.all(adata.obsm[obsm] != adata.obsm[cache_name])

    # use cached
    normalizer(adata, obsm=obsm)
    utils.assert_array_equal(adata.obsm[obsm], adata.obsm[cache_name])
