import unittest
import parameterized
import tempfile
import os
import subprocess

import openproblems
from openproblems.test import utils

utils.ignore_warnings()

TESTDIR = os.path.dirname(os.path.abspath(__file__))


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
    with tempfile.NamedTemporaryFile(suffix=".h5ad") as data_file:
        p = subprocess.run(
            [
                "singularity",
                "exec",
                "docker://singlecellopenproblems/{}".format(method.metadata["image"]),
                "/bin/bash",
                os.path.join(TESTDIR, "singularity_run.sh"),
                TESTDIR,
                "run_test_method.py",
                task_name,
                method.__name__,
                dataset.__name__,
                data_file.name,
            ],
            stderr=subprocess.PIPE,
        )
        assert p.returncode == 0, "Return code {}\n\n{}".format(
            p.returncode, p.stderr.decode("utf-8")
        )
        for metric in task.METRICS:
            image = "docker://singlecellopenproblems/{}".format(
                metric.metadata["image"]
            )
            p = subprocess.run(
                [
                    "singularity",
                    "exec",
                    "docker://singlecellopenproblems/{}".format(
                        metric.metadata["image"]
                    ),
                    "/bin/bash",
                    os.path.join(TESTDIR, "singularity_run.sh"),
                    TESTDIR,
                    "run_test_metric.py",
                    task_name,
                    metric.__name__,
                    data_file.name,
                ],
                stderr=subprocess.PIPE,
            )
            assert p.returncode == 0, "Return code {}\n\n{}".format(
                p.returncode, p.stderr.decode("utf-8")
            )


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
