from . import run

import functools
import openproblems
import os


def git_file_age(filename, branch=None):
    """Get the age of a file's last git commit."""
    git_age = (
        run.run(
            ["git", "log", "-1", '--format="%ad"', "--date=unix", "--", filename],
            return_stdout=True,
        )
        .strip()
        .replace('"', "")
    )
    if git_age == "":
        return 0
    else:
        return int(git_age)


def git_has_diff(filename):
    if isinstance(filename, str):
        filename = [filename]
    diffstat = run.run(
        ["git", "diff", "origin/main", "--shortstat", "--"] + filename,
        return_stdout=True,
    )
    return len(diffstat) > 0


def task_modified(task):
    task_dir = os.path.dirname(task.__file__)
    return git_has_diff(task_dir)


def core_modified():
    return git_has_diff([".", ":^openproblems/tasks"])


@functools.lru_cache(None)
def list_modified_tasks():
    if core_modified():
        return openproblems.TASKS

    return [task for task in openproblems.TASKS if task_modified(task)]
