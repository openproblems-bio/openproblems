import parameterized
import tempfile
import os
import functools
import json
import datetime
import time

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


def docker_image_age(image):
    utils.run(
        ["docker", "images"],
        error_raises=RuntimeError,
        format_error=lambda p: "The Dockerfile for image {} is newer than the "
        "latest push, but Docker is not available. "
        "Return code {}\n\n{}".format(image, p.returncode, p.stdout.decode("utf-8")),
    )
    image_info, returncode = utils.run(
        ["docker", "inspect", "singlecellopenproblems/{}:latest".format(image)],
        return_stdout=True,
        return_code=True,
    )
    if not returncode == 0:
        # image not available
        return 0
    else:
        image_dict = json.loads(image_info)[0]
        created_time = image_dict["Created"].split(".")[0]
        created_timestamp = datetime.datetime.strptime(
            created_time, "%Y-%m-%dT%H:%M:%S"
        ).timestamp()
        return created_timestamp


def docker_push_age(filename):
    try:
        with open(filename, "r") as handle:
            return float(handle.read().strip())
    except FileNotFoundError:
        return 0


@functools.lru_cache(maxsize=None)
def image_requires_docker(image):
    docker_push, dockerfile, requirements = docker_paths(image)
    push_timestamp = docker_push_age(docker_push)
    image_age = utils.git_file_age(dockerfile)
    for req in requirements:
        req_age = utils.git_file_age(req)
        image_age = max(image_age, req_age)
    if push_timestamp > image_age:
        return False
    else:
        if docker_image_age(image) < image_age:
            build_docker(image)
        return True


@functools.lru_cache(maxsize=None)
def cache_singularity_image(image):
    docker_push, _, _ = docker_paths(image)
    push_timestamp = docker_push_age(docker_push)
    image_filename = "{}.sif".format(image)
    image_path = os.path.join(CACHEDIR, image_filename)
    image_age_filename = os.path.join(CACHEDIR, "{}.age.txt".format(image))
    image_age = docker_push_age(image_age_filename)
    if push_timestamp > image_age and os.path.isfile(image_path):
        os.remove(image_path)
    if not os.path.isfile(image_path):
        utils.run(
            [
                "singularity",
                "--verbose",
                "pull",
                "--name",
                image_filename,
                "docker://singlecellopenproblems/{}".format(image),
            ]
        )
        with open(image_age_filename, "w") as handle:
            handle.write(str(time.time()))
    return image_path


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
            "--mount",
            "type=bind,source=/tmp,target=/tmp",
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
