import functools
import git
import openproblems
import pathlib

TASK_MIN_DATASETS = 1
TASK_MIN_METHODS = 3
TASK_MIN_METRICS = 1


@functools.lru_cache()
def get_sha():
    repo = git.Repo(pathlib.Path(openproblems.__path__[0]).parent)
    assert not repo.bare
    return repo.head.commit.hexsha


def task_is_incomplete(task):
    if len(task.DATASETS) < TASK_MIN_DATASETS:
        return True
    non_baseline_methods = [
        method for method in task.METHODS if not method.metadata["is_baseline"]
    ]
    if len(non_baseline_methods) < TASK_MIN_METHODS:
        return True
    if len(task.METRICS) < TASK_MIN_METRICS:
        return True
    return False
