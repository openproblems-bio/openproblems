from .version import __version__

import decorator
import packaging.version
import warnings


@decorator.decorator
def temporary(func, version=None, *args, **kwargs):
    """Decorate a function as a temporary fix.

    Parameters
    ----------
    version : str
        Version after which this function should raise a RuntimeError
    """
    if version is None:
        raise TypeError("temporary() missing 1 required keyword argument: 'version'")
    if packaging.version.parse(__version__) >= packaging.version.parse(version):
        raise RuntimeError(
            "Temporary function {}.{} is temporary and should not be used after version"
            " {} (current version: {})".format(
                func.__module__, func.__name__, version, __version__
            )
        )
    return func(*args, **kwargs)


def future_warning(msg, error_version, error_category, warning_category=FutureWarning):
    """Raise a warning until a specific version, then raise an error."""

    current_version = packaging.version.parse(__version__)
    if current_version < packaging.version.parse(error_version):
        msg += " This will raise a {} in openproblems v{}".format(
            error_category.__name__, error_version
        )
        warnings.warn(msg, warning_category)
    else:
        raise error_category(msg)


def get_members(module):
    """Get all public members from a module."""
    namespace = [attr for attr in dir(module) if not attr.startswith("_")]
    return [getattr(module, attr) for attr in namespace]


def get_callable_members(module):
    """Get all callable public members from a module."""
    return [member for member in get_members(module) if callable(member)]


def get_member_id(member):
    """Get the submodule or function name for a task, dataset, method or metric"""
    return member.__name__.split(".")[-1]
