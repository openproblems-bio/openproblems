# Docker images

By default, all methods and metrics run in the `openproblems` docker image. If you require additional dependencies, you can either add them to an existing docker image, or if this is not possible due to conflicts, add a new one.

To define which image is to be used in a method or metric, simply set the `image` parameter in the method decorator to match the name of the folder containing the Dockerfile (e.g., `image="openproblems-r-base"`).

## Available images

### openproblems

Our base image. Do not add dependencies unless you know what you are doing.

### openproblems-r-base

Our base R image. Do not add dependencies unless you know what you are doing.

### openproblems-r-extras

Our R image that accepts additional dependencies.

To add R packages (CRAN or Bioc), add them to `r_requirements.txt`. Syntax is dictated by [`renv`](https://rstudio.github.io/renv/reference/install.html#examples).

To add Python packages (PyPi or Github), add them to `requirements.txt`. Syntax is dictated by [`pip`](https://packaging.python.org/tutorials/installing-packages/).

### openproblems-python-extras

To add Python packages (PyPi or Github), add them to `requirements.txt`. Syntax is dictated by [`pip`](https://packaging.python.org/tutorials/installing-packages/).

## Adding new images

To add a new image, create a new folder containing the following files:

* `Dockerfile`
* `README.md`
* `requirements.txt` (optional)
* `r_requirements.txt` (optional)

The easiest way to do this is to copy the `openproblems-python-extras` or `openproblems-r-extras` folder.

## Building Docker images locally

If you have Docker installed, you can build containers locally for prototyping. For example, to install the `openproblems` base container, you can run the following.
```
docker build -f docker/openproblems/Dockerfile -t singlecellopenproblems/openproblems .
```
or to update all available Docker images:
```
cd workflow && snakemake -j 10 docker
```
or if you wish to override the automatic change detection,
```
cd workflow && snakemake -j 10 docker_build
```
