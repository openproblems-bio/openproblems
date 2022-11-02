# Guide to Docker images

This guide provides instructions on editing the Docker images used to run `methods`,
`metrics`, and load `datasets` for the Open Problems benchmarking infrastructure.

Note, all images must comply to the [AWS SageMaker Custom Image
Specifications](https://docs.aws.amazon.com/sagemaker/latest/dg/studio-byoi-specs.html).

### Table of Contents  <!-- markdownlint-disable-line MD001 -->

- [About Docker images](#about-docker-images)
- [Available images](#available-images)
  - [openproblems](#openproblems)
  - [openproblems-r-base](#openproblems-r-base)
  - [openproblems-r-extras](#openproblems-r-extras)
  - [openproblems-python-extras](#openproblems-python-extras)
- [Adding a package to the available images](#adding-a-package-to-the-available-images)
- [Adding new images](#adding-new-images)
- [Building Docker images locally](#building-docker-images-locally)
- [Building Docker images through GitHub Actions workflows
  ](#building-docker-images-through-github-actions-workflows)
- [Pulling images from the ECR to your local machine
  ](#pulling-images-from-the-ecr-to-your-local-machine)
- [Running Docker images locally](#running-docker-images-locally)

<!-- Table of contents generated with
[markdown-toc](http://ecotrust-canada.github.io/markdown-toc/) -->

#### Additional resources

- [Dockerfile Reference](https://docs.docker.com/engine/reference/builder/)
  - Documentation from Docker on how to write Dockerfiles
- [SageMaker Studio Custom Image Samples
  ](https://github.com/aws-samples/sagemaker-studio-custom-image-samples/)
  - Example images from AWS designed for compatibility with SageMaker

## About Docker images

By default, all methods, metrics, and dataset loaders run in the `openproblems` docker
image. If you require additional dependencies, you can either add them to an existing
docker image, or if this is not possible due to conflicts, add a new one.

To define which image is to be used in a method or metric, simply set the `image`
parameter in the method decorator to match the name of the folder containing the
Dockerfile (e.g., `image="openproblems-r-base"`).

## Available images

### openproblems

Our base image. Do not add dependencies unless you know what you are doing.

### openproblems-r-base

Our base R image. Do not add dependencies unless you know what you are doing.

### openproblems-r-extras

Our R image that accepts additional dependencies.

To add R packages (CRAN or Bioc), add them to `r_requirements.txt`. Syntax is dictated
by [`renv`](https://rstudio.github.io/renv/reference/install.html#examples).

To add Python packages (PyPi or Github), add them to `requirements.txt`. Syntax is
dictated by [`pip`](https://packaging.python.org/tutorials/installing-packages/).

### openproblems-python-extras

To add Python packages (PyPi or Github), add them to `requirements.txt`. Syntax is
dictated by [`pip`](https://packaging.python.org/tutorials/installing-packages/).

## Adding a package to the available images

Most packages should be able to be added in the Open Problems by editing one of the
available images listed above. If there are conflicting dependencies between the package
you would like to add and the packages already in the available images, follow the
[Adding new images](#adding-new-images) steps below.

Assuming there are no conflicting dependencies, you can simply amend the relevant
`requirements.txt` file in the directory for the Docker image you would like to edit.

1. Select a Docker image to edit. If you're adding a Python package, start with
   `openproblems-python-extras`. If you're adding an R package, start with the
   `openproblems-r-extras`.
2. Edit the relevant `requirements.txt` file.
    - Adding an R package:
        - Edit the `r_requirements.txt` file.
        - The syntax to add a package is defined by [renv
          ](https://rstudio.github.io/renv/reference/install.html#examples).
            - Packages from Bioconductor: `bioc::packagename`
            - Packages from CRAN: `packagename@<version-tag>`
            - Packages from Git: `username/packagename`
        - More complex package installation will require editing the `Dockerfile`.
    - Adding a Python package:
        - Edit the `requirements.txt` file.
        - The syntax to add a package is defined by `pip`
            - Packages from PyPI: `packagename==version`
            - Packages from Git: `git+https://github.com/username/repositoryname`
        - More complex package installation will require editing the `Dockerfile`.
3. Add the `packagename` to the `README.md` file in the directory specifying the Docker
   image. This helps keep track of which packages and versions are installed in each
   Docker image.
4. Commit your changes to the Docker image and push to your fork following the
   instructions in the [Contributing Guide
   ](https://github.com/openproblems-bio/openproblems/blob/main/CONTRIBUTING.md).

## Adding new images

To add a new image, create a new folder containing the following files:

- `Dockerfile`
- `README.md`
- `requirements.txt` (optional)
- `r_requirements.txt` (optional)

The easiest way to do this is to copy the `openproblems-python-extras` or
`openproblems-r-extras` folder.

## Building Docker images locally

If you have Docker installed, you can build containers locally for prototyping. For
example, to install the `openproblems` base container, you can run the following.

```shell
docker build -f docker/openproblems/Dockerfile -t singlecellopenproblems/openproblems .
```

or to update all available Docker images, updating only when necessary:

```shell
cd workflow && snakemake -j 10 docker
```

or if you wish to override the automatic change detection,

```shell
cd workflow && snakemake -j 10 docker_build
```

## Building Docker images through GitHub Actions workflows

Docker images are built by the `run_benchmarks` GitHub Actions workflow on both the base
repository and on forks. As long as you have AWS secrets configured properly for your
repository (see our [Contributing
Guide](https://github.com/openproblems-bio/openproblems/blob/main/CONTRIBUTING.md#submitting-new-features)),
these images will be uploaded to Amazon Web Services [Elastic Container
Registry](https://aws.amazon.com/ecr/) (ECR). You can then download the image locally or
attach to AWS SageMaker Studio.

Once your Run Benchmark has completed successfully, you should see a pane in the GitHub
Actions tab of your fork that looks like this:

![Successful Actions run](https://user-images.githubusercontent.com/8322751/112719533-c508e100-8ecf-11eb-91b0-6f99ccee2e3f.png)

If that workflow failed, you should look at the workflow logs to find the error.

You can find your successfully uploaded images on the ECR. To navigate to the ECR,
search the AWS console for "ECR" and click on "Repositories" and then click on
`openproblems`. You should also see a `nextflow` repository that's used for your
benchmarking backend, but you can ignore that for now.

As you can see below, images uploaded to the ECR have Image Tags in the following format
`openproblems:[first 6 characters of username]-[branch name]-[image name]`. For example,
`danielStrobl` recently pushed his `batch-integration` branch containing a
`openproblems-python37-scgen` image. This is converted to an Image Tag
`daniel-batch-integration-openproblems-python37-scgen`.

![ECR example screen](https://user-images.githubusercontent.com/8322751/112719414-43b14e80-8ecf-11eb-8fe2-5588e42c77c5.png)

## Pulling images from the ECR to your local machine

To pull images from the ECR using `docker pull`, first download and setup the
[`amazon-ecr-credential-helper`](https://github.com/awslabs/amazon-ecr-credential-helper)
using the same AWS secrets that you used to set up your fork repository. With that set
up you can use the following command to pull the image:

```shell
docker pull <aws_account_id>.dkr.ecr.us-west-2.amazonaws.com/openproblems:<Image Tag>
```

If you would like to attach this image to AWS SageMaker, you can follow our [SageMaker
and ECR
tutorial.](https://github.com/openproblems-bio/openproblems/blob/main/SAGEMAKER.md)

You can also pull base images from
[DockerHub](https://hub.docker.com/r/singlecellopenproblems/openproblems):

```shell
docker pull singlecellopenproblems/openproblems-python-extras:latest
```

## Running Docker images locally

To run Docker images on your local machine, you must have `docker` installed. Follow the
Docker guide to [Install Docker](https://docs.docker.com/get-docker/).

Once you've either built Docker images locally or pulled them from ECR or the
[singlecellopenproblems
DockerHub](https://hub.docker.com/r/singlecellopenproblems/openproblems), you can see
installed images using `docker images`.

<!-- markdownlint-disable MD013 -->
```shell
> docker images
REPOSITORY                                                  TAG                                                 IMAGE ID       CREATED        SIZE
singlecellopenproblems/openproblems-python-extras           latest                                              f86e1c5ce9d0   14 hours ago   3.94GB
singlecellopenproblems/openproblems-r-base                  latest                                              f8908c9fb387   21 hours ago   6.36GB
singlecellopenproblems/openproblems-r-extras                latest                                              7e15120bb7ce   5 days ago     4.89GB
singlecellopenproblems/openproblems                         latest                                              14974cbd2f58   5 days ago     2.1GB
490915662541.dkr.ecr.us-west-2.amazonaws.com/openproblems   batch_integration_docker-openproblems               3a1ce37e85f2   6 days ago     2.06GB
```
<!-- markdownlint-enable MD013 -->

You can then run commands within a docker container using `docker run`. Consult the
[Docker documentation](https://docs.docker.com/engine/reference/commandline/run/) to
learn more about the `run` command.

```shell
cd openproblems
docker run \
  -v $(pwd):/usr/src/singlecellopenproblems -v /tmp:/tmp \
  -it singlecellopenproblems/openproblems-python-extras bash
```

You may also specify the docker image by its ID, rather than its name:

```shell
cd openproblems
docker run -v $(pwd):/usr/src/singlecellopenproblems -v /tmp:/tmp -it 90a9110c7d69 bash
```
