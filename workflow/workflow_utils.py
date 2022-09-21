TASK_MIN_DATASETS = 1
TASK_MIN_METHODS = 3
TASK_MIN_METRICS = 1


def task_is_incomplete(task):
    if len(task.DATASETS) < TASK_MIN_DATASETS:
        return True
    if len(task.METHODS) < TASK_MIN_METHODS:
        return True
    if len(task.METRICS) < TASK_MIN_METRICS:
        return True
    return False
