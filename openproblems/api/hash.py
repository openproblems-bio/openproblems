from ..tools.utils import check_version
from . import utils

import hashlib
import importlib
import os
import subprocess

_MODULE = type(os)


def get_module(fun):
    """Get the Python module in which a function is defined."""
    if hasattr(fun, "__wrapped__"):
        fun = fun.__wrapped__
    return fun.__module__


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


def get_context(obj, context=None):
    if context is None:
        context = dict()
    if isinstance(obj, _MODULE):
        module_name = obj.__name__
    else:
        try:
            module_name = get_module(obj)
        except AttributeError:
            # doesn't belong to a module
            return context
    if module_name in context:
        # already done
        pass
    elif not module_name.startswith("openproblems"):
        # use version number of external libs
        version = check_version(module_name.split(".")[0])
        if version != "ModuleNotFound":
            context[module_name] = version
    else:
        module = importlib.import_module(module_name)
        context[module_name] = git_hash(module.__file__)
        for member_name in dir(module):
            member = getattr(module, member_name)
            if callable(member) or isinstance(member, _MODULE):
                context = get_context(member, context=context)
    return context


def hash_dict(context):
    """Convert a SHA256 hash for a Python dictionary of strings."""
    hash = hashlib.sha256()
    for key in sorted(context.keys()):
        hash.update(bytes(key, encoding="utf-8"))
        hash.update(bytes(context[key], encoding="utf-8"))
    return hash.hexdigest()


def get_hash(task_name, function_type, function_name):
    """Get the git commit hash associated with a function."""
    fun = utils.get_function(task_name, function_type, function_name)
    context = get_context(fun)
    hash = hash_dict(context)
    return hash


def main(args):
    """Run the ``hash`` subcommand."""
    hash = get_hash(args.task, args.function_type, args.name)
    return hash
