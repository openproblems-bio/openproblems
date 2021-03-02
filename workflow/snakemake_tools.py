import dateutil.parser
import multiprocessing
import openproblems
import os
import packaging.version
import subprocess

N_THREADS = multiprocessing.cpu_count()
TEMPDIR = ".evaluate"
SCRIPTS_DIR = os.getcwd()
DOCKER_DIR = "/opt/openproblems/scripts/"
RESULTS_DIR = os.path.join(SCRIPTS_DIR, "..", "website", "data", "results")
DOCKER_EXEC = (
    "CONTAINER=$("
    "  docker run -dt --rm"
    '  --mount type=bind,source="{mountdir}",target=/opt/openproblems'
    "  singlecellopenproblems/{{image}}"
    ") bash -c '"
    "  docker exec $CONTAINER /bin/bash /opt/openproblems/scripts/docker_run.sh"
).format(mountdir=os.path.dirname(SCRIPTS_DIR))
try:
    DOCKER_PASSWORD = os.environ["DOCKER_PASSWORD"]
except KeyError:
    DOCKER_PASSWORD = None


def _images(filename):
    docker_dir = os.path.join("..", "docker")
    return [
        os.path.join(docker_dir, image, filename)
        for image in os.listdir(docker_dir)
        if os.path.isdir(os.path.join(docker_dir, image))
    ]


def push_images(wildcards):
    """Get Docker push timestamp for all images."""
    return _images(".docker_push")


def build_images(wildcards):
    """Get Docker build timestamp for all images."""
    return _images(".docker_build")


def pull_images(wildcards):
    """Get Docker pull timestamp for all images."""
    return _images(".docker_pull")


def docker_image_name(wildcards):
    """Get the name of the Docker image required for a task and method/metric."""
    task = getattr(openproblems.tasks, wildcards.task)
    if hasattr(wildcards, "metric"):
        fun = getattr(task.metrics, wildcards.metric)
    elif hasattr(wildcards, "method"):
        fun = getattr(task.methods, wildcards.method)
    else:
        fun = getattr(task.datasets, wildcards.dataset)
    return fun.metadata["image"]


def docker_image_age(image):
    """Get the age of a Docker image."""
    proc = subprocess.run(
        [
            "docker",
            "inspect",
            '-f="{{.Created}}"',
            "singlecellopenproblems/{}".format(image),
        ],
        stdout=subprocess.PIPE,
    )
    date_string = proc.stdout.decode().strip()
    datetime = dateutil.parser.parse(date_string)
    return int(datetime.timestamp())


def docker_file_age(image):
    """Get the age of a Dockerfile."""
    docker_path = "../docker/{}".format(image)
    proc = subprocess.run(
        [
            "git",
            "log",
            "-1",
            '--format="%ad"',
            "--date=unix",
            "--",
            "{}/*".format(docker_path),
        ],
        stdout=subprocess.PIPE,
    )
    return int(proc.stdout.decode().strip())


def version_not_changed():
    """Check that openproblems has not changed version since last build."""
    try:
        with open("../docker/.version", "r") as handle:
            build_version = handle.read().strip()
    except FileNotFoundError:
        return False
    build_version = packaging.version.parse(build_version)
    curr_version = packaging.version.parse(openproblems.__version__)
    if curr_version > build_version:
        return False
    else:
        return True


def docker_image_marker(image):
    """Get the file to be created to ensure Docker image exists from the image name."""
    docker_path = "../docker/{}".format(image)
    docker_push = os.path.join(docker_path, ".docker_push")
    docker_pull = os.path.join(docker_path, ".docker_pull")
    docker_build = os.path.join(docker_path, ".docker_build")
    dockerfile = os.path.join(docker_path, "Dockerfile")
    if version_not_changed() and docker_file_age(docker_push) > docker_image_age(dockerfile):
        # Dockerfile hasn't been changed since last push, pull it
        return docker_pull
    elif DOCKER_PASSWORD:
        # we have the password, let's push it
        return docker_push
    else:
        # new image and we don't have the password, build locally
        return docker_build


def _docker_requirements(image, include_push=False):
    """Get all files to ensure a Docker image is up to date from the image name."""
    docker_dir = "../docker/{}/".format(image)
    dockerfile = os.path.join(docker_dir, "Dockerfile")
    requirements = [dockerfile]
    requirements.extend(
        [
            os.path.join(docker_dir, f)
            for f in os.listdir(docker_dir)
            if f.endswith("requirements.txt")
        ]
    )
    if include_push:
        requirements.append(docker_image_marker(image))
    with open(dockerfile, "r") as handle:
        base_image = next(handle).replace("FROM ", "")
        if base_image.startswith("singlecellopenproblems"):
            base_image = base_image.split(":")[0].split("/")[1]
            requirements.extend(_docker_requirements(base_image, include_push=True))
    return requirements


def docker_requirements(wildcards):
    """Get all files to ensure a Docker image is up to date from wildcards."""
    return _docker_requirements(wildcards.image)


def docker_push(wildcards):
    """Get the file to be created to ensure Docker image exists from wildcards."""
    return docker_image_marker(docker_image_name(wildcards))


def docker_command(wildcards, output):
    """Get the Docker command to be run given a set of wildcards."""
    image = docker_image_name(wildcards)
    return DOCKER_EXEC.format(image=image)
