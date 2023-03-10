from . import utils

import hashlib
import importlib
import json
import os
import random
import scprep
import subprocess
import warnings

_MODULE = type(os)


def _run(args, **kwargs):
    p = subprocess.run(
        args,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        **kwargs,
    )
    if p.returncode != 0:
        raise RuntimeError(p.stderr.decode())
    else:
        return p.stdout.decode()


def get_module(fun):
    """Get the Python module in which a function is defined."""
    if hasattr(fun, "__wrapped__"):
        fun = fun.__wrapped__
    return fun.__module__


def git_hash(obj):
    """Get the git commit hash associated with the latest change to a file."""
    if isinstance(obj, str) and os.path.isfile(obj):
       # if it's a file, run git log to get the hash
        return _run(
            ["git", "log", "-n", "1", "--pretty=format:%H", "--", obj],
            cwd=os.path.dirname(__file__),
        )
    elif hasattr(obj, "__file__"):
        # if it's a module, get the associated file
        return git_hash(obj.__file__)
    elif callable(obj):
        # if it's a function, get the associated module
        return git_hash(importlib.import_module(get_module(obj)))


def docker_token(image_name):
    output = json.loads(
        _run(
            [
                "curl",
                "--silent",
                (
                    f"https://auth.docker.io/token?scope=repository:{image_name}:"
                    "pull&service=registry.docker.io"
                ),
            ]
        )
    )
    return output["token"]


def docker_labels_from_api(image_name, tag="latest"):
    token = docker_token(image_name)
    output = json.loads(
        _run(
            [
                "curl",
                "--silent",
                "--header",
                f"Authorization: Bearer {token}",
                f"https://registry-1.docker.io/v2/{image_name}/manifests/{tag}",
            ]
        )
    )
    v1_compat_output = json.loads(output["history"][0]["v1Compatibility"])
    return v1_compat_output["config"]["Labels"]


def docker_hash(image_name):
    """Get the docker image hash associated with an image."""
    try:
        try:
            return _run(
                [
                    "docker",
                    "inspect",
                    "-f='{{ index .Config.Labels \"bio.openproblems.hash\"}}'",
                    image_name,
                ]
            )
        except (RuntimeError, FileNotFoundError):  # pragma: nocover
            # docker is unavailable or the image is not locally available; use the API
            return docker_labels_from_api(image_name)["bio.openproblems.hash"]
    except Exception:  # pragma: nocover
        warnings.warn(
            "Failed to access docker or the docker API; docker image hash failed. All"
            f" jobs using {image_name} will not be cached."
        )
        return str(random.getrandbits(256))


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
        if isinstance(obj, scprep.run.RFunction) and hasattr(obj, "__r_file__"):
            if obj.__r_file__ not in context:
                context[obj.__r_file__] = git_hash(obj.__r_file__)
        if hasattr(obj, "metadata") and "image" in obj.metadata:
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
