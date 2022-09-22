import warnings


def ignore_warnings():
    """Ignore irrelevant warnings."""
    # warnings.simplefilter("error")

    try:
        import numba

        warnings.filterwarnings("ignore", category=numba.NumbaWarning)
    except ImportError:
        pass
