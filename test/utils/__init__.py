import logging
import tempfile

logging.getLogger("openproblems").setLevel(logging.DEBUG)

TEMPDIR = tempfile.TemporaryDirectory()
