from . import run

import functools
import openproblems
import os


def git_file_age(filename):
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


@functools.lru_cache(None)
def git_setup_remote_base():
    """Set up the `base` remote to point to the main repo."""
    remotes = run.run(["git", "remote"], return_stdout=True).strip().split("\n")
    if "base" not in remotes:
        run.run(
            [
                "git",
                "remote",
                "add",
                "base",
                "https://github.com/openproblems-bio/openproblems",
            ],
            return_stdout=True,
        )
    run.run(["git", "fetch", "base", "main"])


def git_has_diff(filename):
    """Returns true if the filename has changed relative to base/main."""
    if isinstance(filename, str):
        filename = [filename]
    diffstat = run.run(
        ["git", "diff", "base/main", "--shortstat", "--"] + filename,
        return_stdout=True,
    )
    return len(diffstat) > 0


def task_dir(task):
    """Get the base directory of a task"""
    return os.path.dirname(task.__file__)


def task_modified(task):
    """Check if the task has changed relative to base/main."""
    return git_has_diff(task_dir(task))


def core_modified():
    """Check if the core repo has changed relative to base/main.

    We exclude all task directories as well as any md files and the website.
    """
    task_exclusions = [f":^{task_dir(task)}" for task in openproblems.TASKS]
    return git_has_diff([".", ":^*.md", ":^website"] + task_exclusions)


@functools.lru_cache(None)
def list_modified_tasks():
    """List tasks for which testing must be run.

    Return all tasks if the core repo has changed,
    otherwise just those that have changed relative to base/main.
    """
    git_setup_remote_base()
    if core_modified():
        return openproblems.TASKS

    return [task for task in openproblems.TASKS if task_modified(task)]
