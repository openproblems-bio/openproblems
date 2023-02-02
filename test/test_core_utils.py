import numpy as np
import openproblems
import openproblems.tools.utils
import packaging.version
import rpy2.robjects as ro
import sklearn


def test_temporary_version_missing():
    """Test temporary decorator behavior with missing version."""

    @openproblems.utils.temporary
    def test_fn():  # pragma: nocover
        pass

    np.testing.assert_raises_regex(
        TypeError, "missing 1 required keyword argument: 'version'", test_fn
    )


def test_temporary_version_future():
    """Test temporary decorator behavior with version from the future."""
    version = packaging.version.parse(openproblems.__version__)
    if version.minor > 0:
        temp_version = "{}.{}".format(version.major, version.minor - 1)
    else:
        temp_version = "{}.{}".format(version.major - 1, 0)

    @openproblems.utils.temporary(version=temp_version)
    def test_fn():  # pragma: nocover
        pass

    np.testing.assert_raises_regex(
        RuntimeError,
        "Temporary function {}.{} is temporary and should not be used "
        "after version {}".format(test_fn.__module__, test_fn.__name__, temp_version),
        test_fn,
    )


def test_future_warning():
    """Test future warning behavior with version from the past."""
    version = packaging.version.parse(openproblems.__version__)
    error_version = "{}.{}".format(version.major, version.minor + 1)

    np.testing.assert_warns(
        FutureWarning,
        openproblems.utils.future_warning,
        msg="this is a future warning",
        error_version=error_version,
        error_category=ValueError,
    )


def test_future_warning_error():
    """Test future warning behavior with version from the future."""
    version = packaging.version.parse(openproblems.__version__)
    if version.minor > 0:
        error_version = "{}.{}".format(version.major, version.minor - 1)
    else:
        error_version = "{}.{}".format(version.major - 1, 0)

    np.testing.assert_raises_regex(
        ValueError,
        "this is a future warning",
        openproblems.utils.future_warning,
        msg="this is a future warning",
        error_version=error_version,
        error_category=ValueError,
    )


def test_package_version():
    assert openproblems.tools.utils.check_version("scikit-learn") == sklearn.__version__


def test_package_version_missing():
    assert openproblems.tools.utils.check_version("not_a_module") == "ModuleNotFound"


def test_r_package_version():
    version = ro.r("as.character(packageVersion('BiocManager'))")[0]
    assert openproblems.tools.utils.check_r_version("BiocManager") == version


def test_r_package_version_missing():
    assert openproblems.tools.utils.check_r_version("not_a_module") == "ModuleNotFound"
