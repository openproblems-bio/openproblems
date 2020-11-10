import tempfile
import os
import shutil
import atexit


def _make_tempdir():
    tempdir = os.path.join(tempfile.gettempdir(), "openproblems_cache")
    try:
        os.mkdir(tempdir)
    except OSError:
        pass
    return tempdir


TEMPDIR = _make_tempdir()


def _cleanup():
    if os.path.isdir(TEMPDIR):
        try:
            shutil.rmtree(TEMPDIR)
        except (FileNotFoundError, PermissionError):
            pass


atexit.register(_cleanup)


def no_cleanup():
    """Don't delete temporary data files on exit."""
    atexit.unregister(_cleanup)
