from . import asserts, data, docker, name, run, warnings
import logging
import tempfile

logging.getLogger("openproblems").setLevel(logging.DEBUG)

TEMPDIR = tempfile.TemporaryDirectory()
