from . import patch
from . import tasks
from . import tools
from . import utils
from .version import __version__

import logging as _logging

_logging.basicConfig()
log = _logging.getLogger("openproblems")

TASKS = utils.get_members(tasks)
