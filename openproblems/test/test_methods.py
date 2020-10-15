import unittest
import parameterized
import tempfile
import os
import subprocess
import functools

import openproblems
from openproblems.test import utils

utils.ignore_warnings()

TESTDIR = os.path.dirname(os.path.abspath(__file__))
CACHEDIR = os.path.join(tempfile.gettempdir(), ".singularity")
os.environ["SINGULARITY_CACHEDIR"] = CACHEDIR
os.environ["SINGULARITY_PULLFOLDER"] = CACHEDIR


@functools.lru_cache(maxsize=None)
def cache_image(image):
    filename = "{}.sif".format(image)
    p = subprocess.run(
        [
            "singularity",
            "--verbose",
            "pull",
            "--name",
            filename,
            "docker://singlecellopenproblems/{}".format(image),
        ],
        stderr=subprocess.PIPE,
    )
    assert p.returncode == 0, "Return code {}\n\n{}".format(
        p.returncode, p.stderr.decode("utf-8")
    )
    return filename


def singularity_command(image, script, *args):
    return [
        "singularity",
        "--verbose",
        "exec",
        cache_image(image),
        "/bin/bash",
        os.path.join(TESTDIR, "singularity_run.sh"),
        TESTDIR,
        script,
    ] + list(args)


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
            singularity_command(
                method.metadata["image"],
                "run_test_method.py",
                task_name,
                method.__name__,
                dataset.__name__,
                data_file.name,
            ),
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
                singularity_command(
                    metric.metadata["image"],
                    "run_test_metric.py",
                    task_name,
                    metric.__name__,
                    data_file.name,
                ),
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
