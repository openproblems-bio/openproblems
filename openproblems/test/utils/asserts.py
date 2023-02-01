from openproblems.tools.utils import assert_finite  # noqa

import functools
import numpy as np
import pathlib
import scipy.sparse

_REQUEST_HEADERS = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:71.0) "
    "Gecko/20100101 Firefox/71.0"
}
FILEPATH = pathlib.Path(__file__)

_MISSING_DOIS = ["vandermaaten2008visualizing", "hosmer2013applied"]


def assert_array_equal(X, Y):
    """Assert two arrays to be equal, whether sparse or dense."""
    assert X.shape == Y.shape
    if scipy.sparse.issparse(X) and scipy.sparse.issparse(Y):
        X = X.tocsr()
        Y = Y.tocsr()
        X.sort_indices()
        Y.sort_indices()
        np.testing.assert_array_equal(X.data, Y.data)
        np.testing.assert_array_equal(X.indices, Y.indices)
        np.testing.assert_array_equal(X.indptr, Y.indptr)
    else:
        X = np.asarray(X)
        Y = np.asarray(Y)
        np.testing.assert_array_equal(X, Y)


def _response_ok(response):
    if response.ok:
        return True
    if response.status_code == 429:
        # rejected; too many requests
        return True
    return False


@functools.lru_cache(None)
def assert_url_accessible(url):
    import requests

    with requests.head(url, headers=_REQUEST_HEADERS) as response:
        assert _response_ok(response), (url, response.status_code)
    return True


@functools.lru_cache(None)
def _load_bibliography():
    import bibtexparser

    bib_path = FILEPATH.parents[3].joinpath("main.bib")
    with open(bib_path, "r") as handle:
        return bibtexparser.load(handle)


def assert_valid_reference(ref):
    bib = _load_bibliography()
    assert ref in bib.entries_dict
    bibentry = bib.entries_dict[ref]
    if not bibentry["ENTRYTYPE"] == "misc" or ref in _MISSING_DOIS:
        assert "doi" in bibentry
        assert assert_url_accessible(f"https://doi.org/{bibentry['doi']}")
    return True
