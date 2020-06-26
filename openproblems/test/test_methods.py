import unittest
import parameterized
import numbers

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
    adata = dataset(test=True)
    method(adata)
    assert task.checks.check_method(adata)
    for metric in task.METRICS:
        m = metric(adata)
        assert isinstance(m, numbers.Number)
