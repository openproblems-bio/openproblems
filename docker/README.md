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

## Building Docker images through GitHub Actions workflows

Docker images are built by the `run_benchmarks` GitHub Actions workflow on both the base repository and on forks. As long as you have AWS secrets configured properly for your repository (see our [Contributing Guide](https://github.com/singlecellopenproblems/SingleCellOpenProblems/blob/master/CONTRIBUTING.md#submitting-new-features)), these images will be uploaded to Amazon Web Services [Elastic Container Registry](https://aws.amazon.com/ecr/) (ECR). You can then download the image locally or attach to AWS SageMaker Studio.

Once your Run Benchmark has completed successfully, you should see a pane in the GitHub Actions tab of your fork that looks like this:

<img width="800" alt="image" src="https://user-images.githubusercontent.com/8322751/112719533-c508e100-8ecf-11eb-91b0-6f99ccee2e3f.png">

If that workflow failed, you should look at the workflow logs to find the error.

You can find your successfully uploaded images on the ECR. To navigate to the ECR, search the AWS console for "ECR" and click on "Repositories" and then click on `openproblems`. You should also see a `nextflow` repository that's used for your benchmarking backend, but you can ignore that for now.

As you can see below, images uploaded to the ECR have Image Tags in the following format `openproblems:[first 6 characters of username]-[branch name]-[image name]`. For example, `danielStrobel` recently pushed his `batch-integration` branch containing a `openproblems-python37-scgen` image. This is converted to an Image Tag `daniel-batch-integration-openproblems-python37-scgen`.


<img width="800" alt="Untitled" src="https://user-images.githubusercontent.com/8322751/112719414-43b14e80-8ecf-11eb-8fe2-5588e42c77c5.png">

To pull images from the ECR using `docker pull`, first download and setup the [`amazon-ecr-credential-helper`](https://github.com/awslabs/amazon-ecr-credential-helper) using the same AWS secrets that you used to set up your fork repository. With that set up you can use the following command to pull the image:

```
docker pull <aws_account_id>.dkr.ecr.us-west-2.amazonaws.com/openproblems:<Image Tag>
```

If you would like to attach this image to AWS SageMaker, you can follow our [SageMaker and ECR tutorial.](https://github.com/singlecellopenproblems/SingleCellOpenProblems/blob/master/SAGEMAKER.md)
