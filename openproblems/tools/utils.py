import gc
import pkg_resources


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
