import packaging.version
import warnings


def future_warning(msg, error_version, error_category, warning_category=FutureWarning):
    """Raise a warning until a specific version, then raise an error."""
    import openproblems

    current_version = packaging.version.parse(openproblems.__version__)
    if current_version < packaging.version.parse(error_version):
        msg += " This will raise a {} in openproblems v{}".format(
            error_category.__name__, error_version
        )
        warnings.warn(msg, warning_category)
    else:
        raise error_category(msg)


def ignore_warnings():
    """Ignore irrelevant warnings."""
    warnings.filterwarnings(
        "ignore",
        category=FutureWarning,
        message="is_categorical is deprecated and will be removed in a future version."
        "  Use is_categorical_dtype instead",
    )

    try:
        import numba
    except ImportError:
        return

    warnings.filterwarnings("ignore", category=numba.NumbaWarning)
