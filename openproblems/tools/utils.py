import gc
import numpy as np
import pkg_resources
import scipy.sparse


def check_version(pkg):
    """Get the version of a package that may or may not be installed."""
    try:
        return pkg_resources.get_distribution(pkg).version
    except pkg_resources.DistributionNotFound:
        return "ModuleNotFound"


def garbage_collection():
    """Run memory garbage collector.

    Runs gc.collect multiple times to free memory to OS rather than just to Python.
    """
    gc.collect()
    gc.collect()
    gc.collect()


def assert_finite(X):
    """Assert numpy or scipy matrix is finite."""
    X = X.data if scipy.sparse.issparse(X) else X
    assert np.all(np.isfinite(X))
    return True
