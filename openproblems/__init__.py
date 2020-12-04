import logging as _logging
from . import tasks, utils, tools, patch
from .version import __version__

log = _logging.getLogger("openproblems")
TASKS = utils.get_members(tasks)

patch.patch_anndata()
