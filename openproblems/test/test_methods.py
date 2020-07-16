import unittest
import parameterized
import numbers

import openproblems
from openproblems.test import utils

utils.ignore_numba_warnings()


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
    output = method(adata)
    assert output is None
    assert task.checks.check_method(adata)
    for metric in task.METRICS:
        m = metric(adata)
        assert isinstance(m, numbers.Number)
