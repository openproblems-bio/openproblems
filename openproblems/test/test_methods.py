import unittest
import parameterized
import tempfile
import os
import subprocess
import functools
import json
import datetime

import openproblems
from openproblems.test import utils

utils.ignore_warnings()

TESTDIR = os.path.dirname(os.path.abspath(__file__))
BASEDIR = os.path.dirname(os.path.dirname(TESTDIR))
CACHEDIR = os.path.join(os.environ["HOME"], ".singularity")
os.environ["SINGULARITY_CACHEDIR"] = CACHEDIR
os.environ["SINGULARITY_PULLFOLDER"] = CACHEDIR


@functools.lru_cache(maxsize=None)
def image_requires_docker(image):
    docker_path = os.path.join(BASEDIR, "docker", image)
    docker_push = os.path.join(docker_path, ".docker_push")
    dockerfile = os.path.join(docker_path, "Dockerfile")
    if os.path.getmtime(docker_push) > os.path.getmtime(dockerfile):
        return False
    else:
        p = subprocess.run(
            ["docker", "images"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        if not p.returncode == 0:
            raise RuntimeError(
                "The Dockerfile for image {} is newer than the "
                "latest push, but Docker is not available. "
                "Return code {}\n\n{}".format(
                    image, p.returncode, p.stdout.decode("utf-8")
                )
            )
        p = subprocess.run(
            ["docker", "inspect", "singlecellopenproblems/{}:latest".format(image)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if not p.returncode == 0:
            raise RuntimeError(
                "The Dockerfile for image singlecellopenproblems/{image} is newer than the "
                "latest push, but the image has not been built. Build it with "
                "`docker build -f {dockerfile} -t singlecellopenproblems/{image} {basedir}`."
                "Return code {code}\n\n{stderr}".format(
                    image=image,
                    code=p.returncode,
                    stderr=p.stderr.decode("utf-8"),
                    dockerfile=dockerfile,
                    basedir=basedir,
                )
            )
        else:
            image_info = json.loads(p.stdout.decode("utf-8"))[0]
            created_time = image_info["Created"].split(".")[0]
            created_timestamp = datetime.datetime.strptime(
                created_time, "%Y-%m-%dT%H:%M:%S"
            ).timestamp()
            if not created_timestamp > os.path.getmtime(dockerfile):
                raise RuntimeError(
                    "The Dockerfile for image singlecellopenproblems/{image} is"
                    " newer than the latest build. Build it with "
                    "`docker build -f {dockerfile} -t singlecellopenproblems/{image} {basedir}`.".format(
                        image=image, dockerfile=dockerfile, basedir=basedir
                    )
                )
        return True


@functools.lru_cache(maxsize=None)
def cache_singularity_image(image):
    filename = "{}.sif".format(image)
    if not os.path.isfile(os.path.join(CACHEDIR, filename)):
        p = subprocess.run(
            [
                "singularity",
                "--verbose",
                "pull",
                "--name",
                filename,
                "docker://singlecellopenproblems/{}".format(image),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        assert p.returncode == 0, "Return code {}\n\n{}".format(
            p.returncode, p.stdout.decode("utf-8")
        )
    return os.path.join(CACHEDIR, filename)


def singularity_command(image, script, *args):
    return [
        "singularity",
        "--verbose",
        "exec",
        "-B",
        "{}:/opt/openproblems".format(BASEDIR),
        cache_singularity_image(image),
        "/bin/bash",
        "/opt/openproblems/openproblems/test/singularity_run.sh",
        "/opt/openproblems/openproblems/test",
        script,
    ] + list(args)


def cache_docker_image(image):
    p = subprocess.run(
        [
            "docker",
            "run",
            "-dt",
            "--rm",
            "--mount",
            "type=bind,source={},target=/opt/openproblems".format(BASEDIR),
            "singlecellopenproblems/{}".format(image),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert p.returncode == 0, "Return code {}\n\n{}".format(
        p.returncode, p.stderr.decode("utf-8")
    )
    return p.stdout.decode("utf-8")[:12]


def docker_command(image, script, *args):
    container = cache_docker_image(image)
    run_command = [
        "docker",
        "exec",
        container,
        "/bin/bash",
        "/opt/openproblems/openproblems/test/singularity_run.sh",
        "/opt/openproblems/openproblems/test/",
        script,
    ] + list(args)
    stop_command = ["docker", "stop", container]
    return run_command, stop_command


def run_image(image, *args):
    if image_requires_docker(image):
        command, stop_command = docker_command(image, *args)
    else:
        command = singularity_command(image, *args)
    p = subprocess.run(
        command,
        stderr=subprocess.STDOUT,
        stdout=subprocess.PIPE,
    )
    if image_requires_docker(image):
        subprocess.run(stop_command)
    assert p.returncode == 0, "Return code {}\n\n{}".format(
        p.returncode, p.stdout.decode("utf-8")
    )


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
    with tempfile.NamedTemporaryFile(suffix=".h5ad") as data_file:
        run_image(
            method.metadata["image"],
            "run_test_method.py",
            task_name,
            method.__name__,
            dataset.__name__,
            data_file.name,
        )
        for metric in task.METRICS:
            run_image(
                metric.metadata["image"],
                "run_test_metric.py",
                task_name,
                metric.__name__,
                data_file.name,
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
