from openproblems.tools.utils import assert_finite  # noqa

import numpy as np
import scipy.sparse

_REQUEST_HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:71.0) "
    "Gecko/20100101 Firefox/71.0"
}


def assert_array_equal(X, Y):
    """Assert two arrays to be equal, whether sparse or dense."""
    assert X.shape == Y.shape
    if scipy.sparse.issparse(X) and scipy.sparse.issparse(Y):
        X = X.tocsr()
        Y = Y.tocsr()
        np.testing.assert_array_equal(X.data, Y.data)
        np.testing.assert_array_equal(X.indices, Y.indices)
        np.testing.assert_array_equal(X.indptr, Y.indptr)
    else:
        X = np.asarray(X)
        Y = np.asarray(Y)
        np.testing.assert_array_equal(X, Y)


def assert_url_accessible(url):
    import requests

    with requests.head(url, headers=_REQUEST_HEADERS) as response:
        assert response.ok, (url, response.status_code)
    return True
