from . import asserts, cache, data, docker, name, run, warnings
import logging
import tempfile

logging.getLogger("openproblems").setLevel(logging.DEBUG)

TEMPDIR = tempfile.TemporaryDirectory()
