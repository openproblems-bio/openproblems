import parameterized
import tempfile
import os
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


def docker_paths(image):
    docker_path = os.path.join(BASEDIR, "docker", image)
    docker_push = os.path.join(docker_path, ".docker_push")
    dockerfile = os.path.join(docker_path, "Dockerfile")
    requirements = [
        os.path.join(docker_path, f)
        for f in os.listdir(docker_path)
        if f.endswith("requirements.txt")
    ]
    return docker_push, dockerfile, requirements


def build_docker(image):
    _, dockerfile, _ = docker_paths(image)
    utils.run(
        [
            "docker",
            "build",
            "-f",
            dockerfile,
            "-t",
            "singlecellopenproblems/{}".format(image),
            BASEDIR,
        ],
        print_stdout=True,
    )


@functools.lru_cache(maxsize=None)
def image_requires_docker(image):
    docker_push, dockerfile, requirements = docker_paths(image)
    try:
        with open(docker_push, "r") as handle:
            push_timestamp = int(handle.read().strip())
    except FileNotFoundError:
        push_timestamp = 0
    image_age = utils.git_file_age(dockerfile)
    for req in requirements:
        req_age = utils.git_file_age(req)
        image_age = max(image_age, req_age)
    if push_timestamp > image_age:
        return False
    else:
        utils.run(
            ["docker", "images"],
            error_raises=RuntimeError,
            format_error=lambda p: "The Dockerfile for image {} is newer than the "
            "latest push, but Docker is not available. "
            "Return code {}\n\n{}".format(
                image, p.returncode, p.stdout.decode("utf-8")
            ),
        )
        image_info, returncode = utils.run(
            ["docker", "inspect", "singlecellopenproblems/{}:latest".format(image)],
            return_stdout=True,
            return_code=True,
        )
        if not returncode == 0:
            build_docker(image)
        else:
            image_dict = json.loads(image_info)[0]
            created_time = image_dict["Created"].split(".")[0]
            created_timestamp = datetime.datetime.strptime(
                created_time, "%Y-%m-%dT%H:%M:%S"
            ).timestamp()
            if not created_timestamp > os.path.getmtime(dockerfile):
                build_docker(image)
        return True


@functools.lru_cache(maxsize=None)
def cache_singularity_image(image):
    filename = "{}.sif".format(image)
    if not os.path.isfile(os.path.join(CACHEDIR, filename)):
        utils.run(
            [
                "singularity",
                "--verbose",
                "pull",
                "--name",
                filename,
                "docker://singlecellopenproblems/{}".format(image),
            ]
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
    hash = utils.run(
        [
            "docker",
            "run",
            "-dt",
            "--rm",
            "--mount",
            "type=bind,source={},target=/opt/openproblems".format(BASEDIR),
            "singlecellopenproblems/{}".format(image),
        ],
        return_stdout=True,
    )
    return hash[:12]


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
    utils.run(command)
    if image_requires_docker(image):
        utils.run(stop_command)


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
