import unittest
import parameterized
import tempfile

import openproblems
from openproblems.test import utils

utils.ignore_numba_warnings()


@parameterized.parameterized.expand(
    [
        (task, dataset, method)
        for task in openproblems.TASKS
        for dataset in task.DATASETS
        for method in task.METHODS
    ],
    name_func=utils.name_test,
)
def test_method(task, dataset, method):
    task_name = task.__name__.split(".")[-1]
    image = "docker://singlecellopenproblems/{}".format(method.metadata["image"])
    with tempfile.NamedTemporaryFile() as tempfile:
        code = subprocess.call(
            [
                "bash",
                "singularity_run.sh",
                "run_test_method.py",
                task_name,
                method.__name__,
                openproblems.data.TEMPDIR,
                dataset.__name__,
                tempfile,
            ]
        )
        assert code == 0, code
        for metric in task.METRICS:
            image = "docker://singlecellopenproblems/{}".format(
                metric.metadata["image"]
            )
            code = subprocess.call(
                [
                    "bash",
                    "singularity_run.sh",
                    "run_test_metric.py",
                    task_name,
                    metric.__name__,
                    tempfile,
                ]
            )
            assert code == 0, code


@parameterized.parameterized.expand(
    [(method,) for task in openproblems.TASKS for method in task.METHODS],
    name_func=utils.name_test,
)
def test_method_metadata(method):
    assert hasattr(method, "metadata")
    for attr in [
        "method_name",
        "paper_name",
        "paper_url",
        "paper_year",
        "code_url",
        "code_version",
        "image",
    ]:
        assert attr in method.metadata
