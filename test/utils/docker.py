from . import exceptions
from . import git
from . import run

import atexit
import datetime
import decorator
import functools
import inspect
import json
import os
import tempfile
import time
import warnings

TESTDIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASEDIR = os.path.dirname(TESTDIR)
DOCKER_BASEDIR = BASEDIR.replace("/__w", "/home/runner/work")
CACHEDIR = os.path.join(os.environ["HOME"], ".singularity")
os.environ["SINGULARITY_CACHEDIR"] = CACHEDIR
os.environ["SINGULARITY_PULLFOLDER"] = CACHEDIR


def docker_paths(image):
    """Get relevant paths for a Docker image."""
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
    """Build a Docker image."""
    _, dockerfile, _ = docker_paths(image)
    run.run(
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
def docker_available():
    """Check if Docker can be run."""
    returncode = run.run(["docker", "images"], return_code=True)
    return returncode == 0


def docker_image_age(image):
    """Check when the Docker image was built."""
    assert docker_available()
    image_info, returncode = run.run(
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
    """Check when the Docker image was last pushed to Docker Hub."""
    try:
        with open(filename, "r") as handle:
            return float(handle.read().strip())
    except FileNotFoundError:
        return 0


@functools.lru_cache(maxsize=None)
def image_requires_docker(image):
    """Check if a specific image requires Docker or Singularity.

    If the image has been modified more recently than it was pushed to Docker Hub, then
    it should be run in Docker. Otherwise, we use Singularity.
    """
    docker_push, dockerfile, requirements = docker_paths(image)
    push_timestamp = docker_push_age(docker_push)
    git_file_age = git.git_file_age(dockerfile)
    for req in requirements:
        req_age = git.git_file_age(req)
        git_file_age = max(git_file_age, req_age)
    if push_timestamp > git_file_age:
        return False
    else:
        if not docker_available():
            raise RuntimeError(
                "The Dockerfile for image {} is newer than the "
                "latest push, but Docker is not available."
            )
        if docker_image_age(image) < git_file_age:
            import sys

            print(
                "Building {}:\n"
                "Docker push age: {}\n"
                "Docker image modified: {}\n"
                "Docker image age: {}".format(
                    image, push_timestamp, git_file_age, docker_image_age(image)
                ),
                file=sys.stderr,
            )
            sys.stderr.flush()
            build_docker(image)
        return True


@functools.lru_cache(maxsize=None)
def cache_singularity_image(image):
    """Download a Singularity image from Dockerhub."""
    docker_push, _, _ = docker_paths(image)
    push_timestamp = docker_push_age(docker_push)
    image_filename = "{}.sif".format(image)
    image_path = os.path.join(CACHEDIR, image_filename)
    image_age_filename = os.path.join(CACHEDIR, "{}.age.txt".format(image))
    image_age = docker_push_age(image_age_filename)
    if push_timestamp > image_age and os.path.isfile(image_path):
        os.remove(image_path)
    if not os.path.isfile(image_path):
        run.run(
            [
                "singularity",
                "--verbose",
                "pull",
                "--name",
                image_filename,
                "docker://singlecellopenproblems/{}".format(image),
            ],
        )
        with open(image_age_filename, "w") as handle:
            handle.write(str(time.time()))
    return image_path


def singularity_command(image, script, *args):
    """Get the Singularity command to run a script."""
    return [
        "singularity",
        "--verbose",
        "exec",
        "-B",
        "{0}:{0}".format(BASEDIR),
        cache_singularity_image(image),
        "/bin/bash",
        "{}/test/docker_run.sh".format(BASEDIR),
        "{}/test".format(BASEDIR),
        script,
    ] + list(args)


@functools.lru_cache(maxsize=None)
def cache_docker_image(image):
    """Run a Docker image and get the machine ID."""
    tempdir = tempfile.gettempdir()
    hash = run.run(
        [
            "docker",
            "run",
            "-dt",
            "--rm",
            "--user=root",
            "--mount",
            f"type=bind,source={DOCKER_BASEDIR},target={BASEDIR}",
            "--mount",
            f"type=bind,source={tempdir},target={tempdir}",
            f"singlecellopenproblems/{image}",
        ],
        return_stdout=True,
    )
    container = hash[:12]

    def stop():
        run.run(["docker", "stop", container])

    atexit.register(stop)
    return container


def docker_command(image, script, *args):
    """Get the Docker command to run a script."""
    container = cache_docker_image(image)
    run_command = [
        "docker",
        "exec",
        container,
        "/bin/bash",
        f"{BASEDIR}/test/docker_run.sh",
        f"{BASEDIR}/test/",
        script,
    ] + list(args)
    return run_command


def run_image(image, script, *args, timeout=None, retries=0):
    """Run a Python script in a container."""
    if image_requires_docker(image) or docker_available():
        container_command = docker_command
    else:
        container_command = singularity_command
    command = container_command(image, script, *args)
    while True:
        try:
            return run.run(command, timeout=timeout)
        except Exception as e:
            if retries > 0 and not isinstance(e, exceptions.TimeoutError):
                warnings.warn(
                    'Container failed with {}("{}").'.format(type(e).__name__, str(e)),
                    RuntimeWarning,
                )
                retries -= 1
            else:
                raise


@decorator.decorator
def docker_test(func, image=None, timeout=None, retries=0, *args, **kwargs):
    """Run a test function in Docker.

    The function must take only simple objects as arguments
    (i.e. eval(str(args)) == args) and the final argument must be the Docker image.
    """
    if image is None:
        image = args[-1]
    if not image.startswith("openproblems"):
        warnings.warn("Image {} expected to begin with openproblems.".format(image))
    assert eval(str(args)) == args
    assert eval(str(kwargs)) == kwargs
    with tempfile.TemporaryDirectory() as tempdir:
        f = os.path.join(tempdir, "{}.py".format(func.__name__))
        with open(f, "w") as handle:
            in_func = False
            for line in inspect.getsource(func).split("\n"):
                if in_func or line.startswith("def {}(".format(func.__name__)):
                    in_func = True
                    handle.write(line + "\n")
            handle.write("\n")
            handle.write("if __name__ == '__main__':\n")
            handle.write("    import openproblems\n")
            handle.write("    import sys\n")
            handle.write(
                "    sys.path.append('{}')\n".format(
                    os.path.join(DOCKER_BASEDIR, "test")
                )
            )
            handle.write("    {}(*{}, **{})\n".format(func.__name__, args, kwargs))

        run_image(image, f, timeout=timeout, retries=retries)
