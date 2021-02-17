import logging as _logging

from . import patch
from . import tasks
from . import tools
from . import utils
from .version import __version__

_logging.basicConfig()
log = _logging.getLogger("openproblems")

TASKS = utils.get_members(tasks)

patch.patch_anndata()
