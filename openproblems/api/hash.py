from . import utils

import hashlib
import importlib
import os
import scprep
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


def docker_hash(image_name):
    """Get the docker image hash associated with an image."""
    subprocess.run(["docker", "pull", "-q", image_name])
    p = subprocess.run(
        ["docker", "inspect", "--format='{{index .RepoDigests 0}}'", image_name],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=os.path.dirname(__file__),
    )
    if p.returncode != 0:
        raise RuntimeError(p.stderr.decode())
    else:
        return p.stdout.decode()


def get_context(obj, context=None):
    """Return a dictionary describing the context of an object.

    Parameters
    ----------
    obj : callable or module
    context : dict or None
        If not None, adds to the existing dictionary

    Returns
    -------
    context : dict
        Key is a module, value is a hash (internal modules)
        or package version (external modules)
    """
    if context is None:
        context = dict()
    if isinstance(obj, _MODULE):
        module_name = obj.__name__
    else:
        if isinstance(obj, scprep.run.RFunction):
            if obj.__r_file__ not in context:
                try:
                    context[obj.__r_file__] = git_hash(obj.__r_file__)
                except AttributeError:
                    pass
        if hasattr(obj, "metadata"):
            if "image" in obj.metadata:
                image_name = f"singlecellopenproblems/{obj.metadata['image']}"
                if image_name not in context:
                    context[image_name] = docker_hash(image_name)
        try:
            module_name = get_module(obj)
        except AttributeError:
            # doesn't belong to a module
            return context
    if module_name.startswith("openproblems") and module_name not in context:
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
