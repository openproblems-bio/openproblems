from . import utils

import importlib
import os
import subprocess


def get_module(fun):
    """Get the Python module in which a function is defined."""
    if hasattr(fun, "__wrapped__"):
        fun = fun.__wrapped__
    return importlib.import_module(fun.__module__)


def git_hash(file):
    """Get the git commit hash associated with a file."""
    p = subprocess.run(
        ["git", "log", "-n", "1", "--pretty=format:%H", "--", file],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=os.path.dirname(__file__),
    )
    if p.returncode != 0:
        raise RuntimeError(p.stderr.decode())
    else:
        return p.stdout.decode()


def get_hash(task_name, function_type, function_name):
    """Get the git commit hash associated with a function."""
    fun = utils.get_function(task_name, function_type, function_name)
    module = get_module(fun)
    hash = git_hash(module.__file__)
    return hash


def main(args):
    """Run the ``hash`` subcommand."""
    hash = get_hash(args.task, args.function_type, args.name)
    return hash
