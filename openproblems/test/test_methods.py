import unittest
import parameterized

import openproblems


@parameterized.parameterized(
    [
        (task, dataset, method)
        for task in openproblems.TASKS
        for dataset in task.DATASETS
        for method in task.METHODS
    ]
)
def test_method(task, dataset, method):
    adata = dataset()
    method(adata)
    for metric in task.METRICS:
        m = metric(adata)
