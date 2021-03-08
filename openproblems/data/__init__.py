import atexit
import logging
import os
import shutil
import tempfile

log = logging.getLogger("openproblems")


def _make_tempdir():
    tempdir = os.path.join(tempfile.gettempdir(), "openproblems_cache")
    try:
        os.mkdir(tempdir)
        log.debug("Created data cache directory")
    except OSError:
        log.debug("Data cache directory exists")
    return tempdir


TEMPDIR = _make_tempdir()


def _cleanup():
    if os.path.isdir(TEMPDIR):
        try:
            shutil.rmtree(TEMPDIR)
            log.debug("Removed data cache directory")
        except FileNotFoundError:
            log.debug("Data cache directory does not exist, cannot remove")
        except PermissionError:
            log.debug("Missing permissions to remove data cache directory")


atexit.register(_cleanup)


def no_cleanup():
    """Don't delete temporary data files on exit."""
    atexit.unregister(_cleanup)
