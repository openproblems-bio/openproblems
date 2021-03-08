from .version import __version__

import decorator
import packaging.version


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
            "Temporary function {}.{} is temporary and should not be used "
            "after version {} (current version: {})".format(
                func.__module__, func.__name__, version, __version__
            )
        )
    return func(*args, **kwargs)


def get_members(module):
    """Get all public members from a module."""
    namespace = [attr for attr in dir(module) if not attr.startswith("_")]
    return [getattr(module, attr) for attr in namespace]


def get_callable_members(module):
    """Get all callable public members from a module."""
    return [member for member in get_members(module) if callable(member)]
