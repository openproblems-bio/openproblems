# Docker images

By default, all methods and metrics run in the `openproblems` docker image. If you require additional dependencies, you can either add them to an existing docker image, or if this is not possible due to conflicts, add a new one.

To define which image is to be used in a method or metric, simply set the `image` parameter in the method decorator to match the name of the folder containing the Dockerfile (e.g., `image="openproblems-r-base"`).

## Available images

### openproblems

Our base image. Do not add dependencies unless you know what you are doing.

### openproblems-r-base

Our base R image. Do not add dependencies unless you know what you are doing.

### openproblems-r-extras

Our R image that accepts additional dependencies. To add R packages (CRAN or Bioc), add them to `r_requirements.txt`. To add Python packages (PyPi or Github), add them to `requirements.txt`.

### openproblems-python-extras

Our Python image that accepts additional dependencies. To add Python packages (PyPi or Github), add them to `requirements.txt`.

## Adding new images

To add a new image, create a new folder containing the following files:

* `Dockerfile`
* `requirements.txt`
* `README.md`

The easiest way to do this is to copy the `openproblems-python-extras` folder.
