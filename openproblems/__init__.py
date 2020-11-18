from . import tasks, utils, tools, patch
from .version import __version__

TASKS = utils.get_members(tasks)

patch.patch_anndata()
