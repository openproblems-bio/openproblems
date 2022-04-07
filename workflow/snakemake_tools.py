import datetime
import functools
import multiprocessing
import openproblems
import os
import packaging.version
import subprocess
import sys
import time
import warnings

N_THREADS = multiprocessing.cpu_count()
TEMPDIR = ".evaluate"
SCRIPTS_DIR = os.getcwd()
RESULTS_DIR = os.path.join(SCRIPTS_DIR, "..", "website", "data", "results")
TEST_DIR = os.path.join(SCRIPTS_DIR, "..", "test")
IMAGES_DIR = os.path.join("..", "docker")
VERSION_FILE = os.path.join(IMAGES_DIR, ".version")
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

sys.path.append(TEST_DIR)


def _images(filename):
    return [
        os.path.join(IMAGES_DIR, image, filename)
        for image in os.listdir(IMAGES_DIR)
        if os.path.isdir(os.path.join(IMAGES_DIR, image))
    ]


def image_markers(wildcards):
    """Get the appropriate marker for each image."""
    return [
        docker_image_marker(image)
        for image in os.listdir(IMAGES_DIR)
        if os.path.isdir(os.path.join(IMAGES_DIR, image))
    ]


def push_images(wildcards):
    """Get Docker push timestamp for all images."""
    images = _images(".docker_push")
    return images


def build_images(wildcards):
    """Get Docker build timestamp for all images."""
    return _images(".docker_build")


def pull_images(wildcards):
    """Get Docker pull timestamp for all images."""
    return _images(".docker_pull")


def refresh_images(wildcards):
    """Get Docker refresh timestamp for all images."""
    return _images(".docker_refresh")


def update_images(wildcards):
    """Get Docker update timestamp for all images."""
    return _images(".docker_update")


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


def docker_image_exists(image, local=True):
    """Check if a Docker image exists."""
    if local:
        proc = subprocess.run(
            [
                "docker",
                "inspect",
                "singlecellopenproblems/{}".format(image),
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    else:
        env = os.environ.copy()
        env["DOCKER_CLI_EXPERIMENTAL"] = "enabled"
        proc = subprocess.run(
            [
                "docker",
                "manifest",
                "inspect",
                "singlecellopenproblems/{}".format(image),
            ],
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    return proc.returncode == 0


def _docker_base(image):
    """Get the base image from a Dockerfile."""
    dockerfile = os.path.join(IMAGES_DIR, image, "Dockerfile")
    with open(dockerfile, "r") as handle:
        base_image = next(handle).replace("FROM ", "")
        if base_image.startswith("singlecellopenproblems"):
            base_image = base_image.split(":")[0].split("/")[1]
            return base_image
        else:
            return None


def docker_image_age(image, pull_on_error=True):
    """Get the age of a Docker image."""
    proc = subprocess.run(
        [
            "docker",
            "inspect",
            '-f="{{.Created}}"',
            "singlecellopenproblems/{}".format(image),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    date_string = proc.stdout.decode().strip().replace('"', "").split(".")[0]
    try:
        date_datetime = datetime.datetime.strptime(date_string, "%Y-%m-%dT%H:%M:%S")
        return int(date_datetime.timestamp())
    except ValueError:
        if pull_on_error and docker_image_exists(image, local=False):
            print(
                "Pulling image singlecellopenproblems/{}".format(image), file=sys.stderr
            )
            subprocess.run(
                [
                    "docker",
                    "pull",
                    "--quiet",
                    "singlecellopenproblems/{}".format(image),
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            return docker_image_age(image, pull_on_error=False)
        elif date_string == "":
            warnings.warn(
                "Docker image singlecellopenproblems/{} not found; "
                "assuming needs rebuild. If you think this message is in error, "
                "you can fix this by running `snakemake -j 1 docker_pull`".format(image)
            )
            return -1
        else:
            raise


def git_file_diff(filename):
    """Check if there are (un)staged changes to a file."""
    proc = subprocess.run(
        [
            "git",
            "status",
            "--porcelain",
            "--untracked-files=no",
            filename,
        ],
        stdout=subprocess.PIPE,
    )
    result = proc.stdout.decode().strip()
    return result != ""


def git_file_age(filename):
    """Get the age of a file on git."""
    if git_file_diff(filename):
        # there exist (un)staged changes, modified now
        return int(time.time())
    # check when the last committed changes occurred
    proc = subprocess.run(
        [
            "git",
            "log",
            "-1",
            '--format="%ad"',
            "--date=unix",
            "--",
            filename,
        ],
        stdout=subprocess.PIPE,
    )
    result = proc.stdout.decode().strip().replace('"', "")
    try:
        return int(result)
    except ValueError:
        if result == "":
            warnings.warn(
                "Files {}/{} not found on git; assuming unchanged.".format(
                    os.getcwd(), filename
                )
            )
            return 0
        else:
            raise


def docker_file_age(image):
    """Get the age of a Dockerfile."""
    docker_path = os.path.join(IMAGES_DIR, image)
    result = git_file_age("{}/Dockerfile".format(docker_path))
    # check for changes to base image
    base_image = _docker_base(image)
    if base_image is not None:
        result = max(result, docker_file_age(base_image))
    return result


def version_not_changed():
    """Check that openproblems has not changed version since last build."""
    try:
        with open(VERSION_FILE, "r") as handle:
            build_version = handle.read().strip()
    except FileNotFoundError:
        return False
    build_version = packaging.version.parse(build_version)
    curr_version = packaging.version.parse(openproblems.__version__)
    if curr_version > build_version:
        return False
    else:
        return True


def format_timestamp(ts):
    """Format a unix timestamp as a string."""
    return datetime.datetime.fromtimestamp(ts).isoformat()


def docker_imagespec_changed(image, dockerfile):
    """Check if the Dockerfile has changed

    If working with a locally built image, simply check if the local image
    was built more recently than the dockerfile changed.

    If working with a github actions-built image, check if there is any diff
    between the Dockerfile and base/main
    """
    with subprocess.run(
        [
            "docker",
            "inspect",
            "-f='{{ index .Config.Labels \"bio.openproblems.build\"}}'",
            "singlecellopenproblems/{}".format(image),
        ],
        stdout=subprocess.PIPE,
    ) as proc:
        build_type = proc.stdout.strip()
    if build_type == "github_actions":
        import utils.git

        has_diff = utils.git.git_has_diff(dockerfile)
        msg = "changed with respect" if has_diff else "identical"
        print(
            "{}: Dockerfile {} to base/main".format(image, msg),
            file=sys.stderr,
        )
        return has_diff
    elif build_type == "local":
        dockerfile_timestamp = docker_file_age(image)
        docker_image_timestamp = docker_image_age(image)
        print(
            "{}: Dockerfile changed {}; Docker image refreshed {}".format(
                image,
                format_timestamp(dockerfile_timestamp),
                format_timestamp(docker_image_timestamp),
            ),
            file=sys.stderr,
        )
        return dockerfile_timestamp > docker_image_timestamp
    else:
        raise NotImplementedError(f"Unrecognized build type {build_type}")


@functools.lru_cache(maxsize=None)
def docker_image_marker(image, refresh=True):
    """Get the file to be created to ensure Docker image exists from the image name."""
    docker_path = os.path.join(IMAGES_DIR, image)
    # possible outputs
    docker_refresh = os.path.join(docker_path, ".docker_refresh")
    dockerfile = os.path.join(docker_path, "Dockerfile")
    no_change = docker_refresh if refresh else dockerfile
    no_change_text = "refreshing source code only" if refresh else "no change"
    docker_build = os.path.join(docker_path, ".docker_build")

    # inputs to conditional logic
    local_imagespec_changed = docker_imagespec_changed(image, dockerfile)
    local_codespec_changed = not version_not_changed()
    if local_codespec_changed:
        print(
            "Code version changed: {}".format(openproblems.__version__), file=sys.stderr
        )
    if local_imagespec_changed or local_codespec_changed:
        # spec has changed, let's rebuild
        print("{}: rebuilding".format(image), file=sys.stderr)
        requirement_file = docker_build
    elif docker_image_exists(image, local=True) or docker_image_exists(
        image, local=False
    ):
        # docker exists on dockerhub and no changes required
        print("{}: {}".format(image, no_change_text), file=sys.stderr)
        requirement_file = no_change
    else:
        # image doesn't exist anywhere, need to build it
        print("{}: building".format(image), file=sys.stderr)
        requirement_file = docker_build
    sys.stderr.flush()
    return requirement_file


def _docker_requirements(image, include_self=False, refresh=True):
    """Get all files to ensure a Docker image is up to date from the image name."""
    docker_dir = os.path.join(IMAGES_DIR, image)
    dockerfile = os.path.join(docker_dir, "Dockerfile")
    requirements = [dockerfile]
    requirements.extend(
        [
            os.path.join(docker_dir, f)
            for f in os.listdir(docker_dir)
            if f.endswith("requirements.txt")
        ]
    )
    if include_self:
        requirements.append(docker_image_marker(image, refresh=refresh))
    base_image = _docker_base(image)
    if base_image is not None:
        requirements.extend(
            _docker_requirements(base_image, include_self=True, refresh=refresh)
        )
    return requirements


def docker_requirements(wildcards):
    """Get all files to ensure a Docker image is up to date from wildcards."""
    return _docker_requirements(wildcards.image)


def docker_update_requirements(wildcards):
    """Get all files to ensure a Docker image is built and up to date from wildcards."""
    return _docker_requirements(wildcards.image, include_self=True, refresh=True)


def docker_push_requirements(wildcards):
    """Get all files to ensure a Docker image is ready to push from wildcards."""
    return _docker_requirements(wildcards.image, include_self=True, refresh=False)


def docker_push(wildcards):
    """Get the file to be created to ensure Docker image exists from wildcards."""
    return docker_image_marker(docker_image_name(wildcards))


def docker_command(wildcards, output):
    """Get the Docker command to be run given a set of wildcards."""
    image = docker_image_name(wildcards)
    return DOCKER_EXEC.format(image=image)


for image in _images(""):
    for filename in [
        ".docker_push",
        ".docker_pull",
        ".docker_build",
        ".docker_refresh",
    ]:
        try:
            os.remove(os.path.join(image, filename))
        except FileNotFoundError:
            pass
