from .. import TASKS
from . import utils


def get_tasks():
    """Get the name of the each of the openproblems tasks."""
    tasks = [utils.module_to_str(task) for task in TASKS]
    return tasks


def main(args):
    """Run the ``tasks`` subcommand."""
    tasks = get_tasks()
    return tasks
