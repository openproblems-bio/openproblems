import gc
import numpy as np
import pkg_resources
import rpy2.rinterface_lib.embedded
import scipy.sparse
import scprep.run

_check_r_version = scprep.run.RFunction(
    args="pkg", body="as.character(packageVersion(pkg))"
)


def check_version(pkg):
    """Get the version of a Python package that may or may not be installed."""
    try:
        return pkg_resources.get_distribution(pkg).version
    except pkg_resources.DistributionNotFound:
        return "ModuleNotFound"


def check_r_version(pkg):
    """Get the version of an R package that may or may not be installed."""
    try:
        return _check_r_version(pkg)[0]
    except rpy2.rinterface_lib.embedded.RRuntimeError as e:
        if "there is no package called" in str(e):
            return "ModuleNotFound"
        else:
            raise


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
