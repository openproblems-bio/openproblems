# Guide to GitHub and AWS SageMaker for the Open Problems Jamboree

[![Video](https://img.shields.io/static/v1?label=YouTube&message=Watch%20tutorial&color=red&logo=youtube)](https://www.youtube.com/watch?v=mNu8-KR7UFY)

This document is a guide to using the Open Problems in Single-Cell Analysis GitHub repository to edit Docker containers and prototype code in AWS SageMaker. In this guide, you will perform the following steps

### Steps

1. Fork the SingleCellOpenProblems GitHub repository
2. Create a custom Dockerfile
3. Push this image to the Elastic Container Registry using GitHub Actions workflows
4. Attach this image to an AWS SageMaker Studio domain
5. Launch a notebook using the custom image

In addition to this guide, we've produced a video tutorial, [**Using GitHub and AWS for Open Problems in Single-Cell Analysis**](https://www.youtube.com/watch?v=mNu8-KR7UFY), which walks through every step of the process. In the video, we provide details that are helpful for troubleshooting if you're having issues. The video also performs all of the steps in the process from start to finish so you can see what needs to happen at each step. However, because of the troubleshooting tips and included wait times, the video is 40 minutes long. To make it quicker to get started, we've written out the steps in the process here and included time-stamped links to the video in case you'd like to see a specific step performed live.

You also might find it helpful to consult AWS Documentation for bringing your own image to AWS SageMaker. You have two options using either the command-line interface (CLI) or a web interface (Console):
* [Custom Image with SageMaker Studio (CLI)](https://docs.aws.amazon.com/sagemaker/latest/dg/studio-byoi-create-sdk.html)
* [Custom Image with SageMaker Studio (Console)](https://docs.aws.amazon.com/sagemaker/latest/dg/studio-byoi.html).

There is a 1:1 correspondence between the steps to set up SageMaker using the CLI and the Console UI. We use the Console in our video and this guide, but feel free to use the CLI if you are more comfortable with it.

### Table of Contents
- [Introduction to the GitHub repository](#introduction-to-the-github-repository)
  * [Why Docker containers?](#why-docker-containers)
  * [Prototyping within Docker containers](#prototyping-within-docker-containers)
- [Getting started and forking the GitHub repository](#getting-started-and-forking-the-github-repository)
- [Editing a Dockerfile](#editing-a-dockerfile)
- [Find your Docker container on the Elastic Container Registry](#find-your-docker-container-on-the-elastic-container-registry)
- [Attach your Image to SageMaker Studio](#attach-your-image-to-sagemaker-studio)
- [Add user to SageMaker Studio](#add-user-to-sagemaker-studio)
- [Open SageMaker Studio and Launch a Notebook using a Custom Image](#open-sagemaker-studio-and-launch-a-notebook-using-a-custom-image)
  * [Selecting an instance type](#selecting-an-instance-type)
  * [Failed to start kernel (image does not exist)](#failed-to-start-kernel-image-does-not-exist)

<!-- Table of contents generated with [markdown-toc](http://ecotrust-canada.github.io/markdown-toc/) -->


## Introduction to the GitHub repository
[Watch this section of the tutorial starting at [0:00]](https://www.youtube.com/watch?v=mNu8-KR7UFY&t=0s)

Open Problems in Single-Cell Analysis is a project to aggregate and benchmark solutions to formalized problems in single-cell analysis.

To facilitate comparing methods designed for a particular task, we've created a platform that runs benchmarks on AWS servers. This allows individuals to upload code to a central repository and have their method benchmarked using a set of standardized datasets without needing access to compute resources to run the full pipeline.

### Why Docker containers?

Because methods developers use different packages with their own sets dependencies and language requirements, we run each method, metric, and dataset loader in a Docker container. You can configure these containers to use a specific operating system, programming language, and any dependencies needed for a method. Each Docker containers is specified by a `Dockerfile` within the [`docker/`](https://github.com/singlecellopenproblems/SingleCellOpenProblems/tree/master/docker) directory of the Open Problems GitHub repository. Consult the [Docker documentation](https://docs.docker.com/engine/reference/builder/) for guidance on creating a `Dockerfile`.

You can find instructions for adding a new Docker image in our [Docker README.md](https://github.com/singlecellopenproblems/SingleCellOpenProblems/blob/master/docker/README.md).


### Prototyping within Docker containers
To facilitate prototyping within Docker containers, we've set up the Open Problems workflows to automatically compile containers and upload them to the Amazon Web Services (AWS) [Elastic Container Registry](https://aws.amazon.com/ecr/).

You can attach these Images to an AWS SageMaker Domain. From there, you can launch Jupyter notebooks using your custom image. You can then test out code contributions within the image before committing to the GitHub repository.

## Getting started and forking the GitHub repository
[Watch this section of the tutorial starting at [1:20]](https://www.youtube.com/watch?v=mNu8-KR7UFY&t=80s)

The goal for this section is to create a fork, activate GitHub Actions, configure your AWS secrets, edit the GitHub repository locally, and push to your fork.

Before you get started, please read our [Contributor Guide](https://github.com/singlecellopenproblems/SingleCellOpenProblems/blob/master/CONTRIBUTING.md). It contains important up-to-date information about how to best contribute to Open Problems.

Detailed instructions for the forking process can be found in the [Submitting new features](https://github.com/singlecellopenproblems/SingleCellOpenProblems/blob/master/CONTRIBUTING.md#submitting-new-features) section of our Contributor Guide.

## Editing a Dockerfile
[Watch this section of the tutorial starting at [10:08]](https://youtu.be/mNu8-KR7UFY?t=607)

All the Open Problems Docker containers are specified by Dockerfiles in the `docker/` folder.

To edit your custom Docker image on your fork, follow the instructions for editing Docker images in our [Docker Guide](https://github.com/singlecellopenproblems/SingleCellOpenProblems/blob/master/docker/README.md).

Note, if you're interested in creating a container for R, you need to inherit from the `openproblems-r-base` image.

Additional resources
* [SageMaker Studio Custom Image Samples](https://github.com/aws-samples/sagemaker-studio-custom-image-samples) - This GitHub contains example images that already work with SageMaker Studio
* [Bringing your own R environment to Amazon SageMaker Studio](https://aws.amazon.com/blogs/machine-learning/bringing-your-own-r-environment-to-amazon-sagemaker-studio/) - This AWS ML Blog details the start to finish process of creating an R environment for SageMaker

## Find your Docker container on the Elastic Container Registry
[Watch this section starting at [15:08]](https://www.youtube.com/watch?v=mNu8-KR7UFY&t=905s)

Once the *Run Benchmarks* workflow has successfully completed the Upload Docker job, your docker images are uploaded to the AWS [Elastic Container Registry (ECR)](https://aws.amazon.com/ecr/).

You can navigate to the ECR from the AWS Console. You should have received an email with your AWS User Credentials that includes a Console Login URL. Enter your username and password to login. You can then search for the `Elastic Container Registry` service in the search bar at the top of the screen.

Next, follow the steps to locate your image on the ECR in the [Building Docker Images using Github Actions workflows](https://github.com/singlecellopenproblems/SingleCellOpenProblems/blob/master/docker/README.md#building-docker-images-through-github-actions-workflows) section of our Docker guide.

## Attach your Image to SageMaker Studio
[Watch this section of the tutorial starting at [18:50]](https://www.youtube.com/watch?v=mNu8-KR7UFY&t=1130s)

To use your image in a SageMake Studio notebook, you need to attach the image to the domain running the SageMaker app.

First, navigate to the SageMaker Studio control panel from the AWS Console. You should have received an email with your AWS User Credentials that includes a Console Login URL. You can then search for the `Amazon SageMaker` service in the search bar at the top of the screen.

Next, follow the steps in the AWS Documentation to [Attach a Custom Image to an Existing SageMaker Studio Domain](https://docs.aws.amazon.com/sagemaker/latest/dg/studio-byoi-attach.html).

If you have already launched SageMaker Studio, you will need to restart your SageMaker Studio app so it can see the new image:
1. Follow the steps in the [Shut Down SageMaker Resources](https://docs.aws.amazon.com/sagemaker/latest/dg/notebooks-run-and-manage-shut-down.html) tutorial.
2. Follow the steps below to [Open SageMaker Studio](https://github.com/singlecellopenproblems/SingleCellOpenProblems/blob/master/SAGEMAKER.md#open-sagemaker-studio-and-launch-a-notebook-using-a-custom-image).

A few caveats to watch out for during this step:
* The Image name should be the Image Tag from the ECR. **Use only lower case letters and dashes** (should match RegEx `[a-z\-]`). No underscores or other special characters are allowed.
* When editing the Image name and Image display name, you must click outside the box to validate the text.
* The `IAM role` should be the `AmazonSageMaker-ExecutionRole` available within the dropdown.
* The EFS Mount Path should be `/home/sagemaker-user`. As long as your image inherits from the base `openproblems` image, this folder should exist.
* The Kernel name is set within the Docker Image. If you’re using one of the Python images, your kernel name should be `python3`. If you’re using an R image and want an R kernel, the image should be `ir`. These Kernel names must be available to Jupyter KernelGateway. More information about Jupyter Kernels can be found in the [Making kernels for Jupyter](https://jupyter-client.readthedocs.io/en/stable/kernels.html) documentation.
* If you see an error after you click Submit stating that the domain already being updated, just wait a few seconds. Someone else tried to attach an image to the same domain at the same time as you, and the domain can only process one request at a time.

Additional resources
* [Custom SageMaker image specifications](https://docs.aws.amazon.com/sagemaker/latest/dg/studio-byoi-specs.html) - specifications for SageMaker Studio Images
* [Field Notes: Accelerate Research with Managed Jupyter on Amazon SageMaker](https://aws.amazon.com/blogs/architecture/field-notes-accelerate-research-with-managed-jupyter-on-amazon-sagemaker/) - a start to finish blog on creating a SageMaker environment with custom images
* [SageMaker Custom Images - Issue #9](https://github.com/aws-samples/sagemaker-studio-custom-image-samples/issues/9) - an issue @dburkhardt encountered while creating images for this project. SageMaker does some wonky modifications of the images before launching a notebook. This may be a hangup for you if you try to install a custom Python environment.


## Add user to SageMaker Studio
[Watch this section of the tutorial starting at [27:16]](https://www.youtube.com/watch?v=mNu8-KR7UFY&t=1635s)

You only need to do this step if you haven't already created a user on the SageMaker Studio console for the domain you're using.

To add a new user, go to the SageMaker Studio Control Panel and click Add User in the upper right hand corder.

Your user name should be the first letter of the your name followed by your last name. E.g. Wes Lewis becomes `wlewis`.

Execution role should be `AmazonSageMaker-ExecutionRole`. This is the same as used when Attaching the image.

The user configuration page should look like this:

<img src="https://i.imgur.com/BisOfAm.png" width="600px">


To create the user, hit Submit. If you see an error about the domain status, wait a few seconds and try again. One only user can be added to a domain at a given time.

## Open SageMaker Studio and Launch a Notebook using a Custom Image
[Watch this section of the tutorial starting at [28:06]](https://www.youtube.com/watch?v=mNu8-KR7UFY&t=1685s)


To launch SageMaker, click on the Open Studio button next to your user name.

<img src="https://i.imgur.com/nLgtrmI.png" width="150px">

Note, you can access anyone's Studio to facilitate collaboration. Please be aware that this means you can also overwrite their notebooks. Be careful.

If you haven’t launched the Studio before, you should see a loading screen for SageMaker Studio. It may take 2-3 minutes for the server to load.

<img src="https://i.imgur.com/gH9EFwm.png" width="600px">

Once it’s working, you should see a Jupyter Lab interface. For more information on the UI, please consult the AWS [Guide to the SageMaker Studio Interface](https://docs.aws.amazon.com/sagemaker/latest/dg/studio-ui.html).

To launch a Notebook with your custom image:
1. Scroll down to the "Notebooks and Compute Resources"
2. Click on the Select a SageMaker Image dropdown
3. Select your custom image

<img src="https://i.imgur.com/9RY2Ma1.png" width="600px">

4. Select either a Notebook or a Console depending on which you'd like to use. You can also launch an Image terminal, which will start a Bash environment from the Docker image.
5. If you click Notebook, SageMaker will launch a Jupyter Notebook using your custom image. You can see the Kernel being used in the upper right corner of the notebook.

<img src="https://user-images.githubusercontent.com/8322751/112723972-8121d600-8ee7-11eb-9f99-c28929d546c0.png" width="400px">

7. While the image is loading, you will see "Unknown" next to the Kernel name. This will change once SageMaker has requisitioned a virtual instance for your notebook to run in. After the notebook has spun up, it will display information about the vCPU and RAM available in your instance. You may click on this information to select a new instance type.

You can now start programming in this Notebook using the installed packages!

<img src="https://i.imgur.com/eb6aJBG.png" width="600px">

Note, you can also start up R notebooks by selecting an image using the R ipython kernel.

<img src="https://user-images.githubusercontent.com/8322751/112722393-af031c80-8edf-11eb-9a65-3233ffe88fd7.png" width="600px">

8. When you're done, be sure to [Shut Down SageMaker Resources](https://docs.aws.amazon.com/sagemaker/latest/dg/notebooks-run-and-manage-shut-down.html) to avoid incurring excess costs.

### Selecting an instance type

When you open a new notebook for the first time, you are assigned a default Amazon Elastic Compute Cloud (Amazon EC2) instance type to run the notebook. When you open additional notebooks on the same instance type, the notebooks run on the same instance as the first notebook, even if the notebooks use different kernels.

We've selected three instances to use during the Jamboree. Note, it is possible to select any kind of instance, but please only select from the following instances. If you need access to a different instance type, please contact an `@organizer` on Discord.

|   Instance Type   | vCPU   | Memory    | GPU  |  $/hr     | Intended Use                              |
|-------------------|:------:|:---------:|:----:|:---------:|-------------------------------------------|
|    `ml.t3.medium` |   `2`  |  `4 GiB`  |  `0` | `$0.0582` | General prototyping with test datasets    |
|   `ml.t3.2xlarge` |   `8`  |  `32 GiB` |  `0` | `$0.4659` | Production prototyping with full datasets |
| `ml.g4dn.2xlarge` |   `8`  |  `32 GiB` |  `1` | `$1.0528` | Prototyping requiring a GPU               |

To change your instance, follow the [Change Instance Type](https://docs.aws.amazon.com/sagemaker/latest/dg/notebooks-run-and-manage-switch-instance-type.html) tutorial from AWS.


### Failed to start kernel (image does not exist)

If you see the following error:

<img src="https://user-images.githubusercontent.com/8322751/112722196-65fe9880-8ede-11eb-83bb-72716866b411.png" width="600px">

Don't fret! This means that the image you were using has been deleted from ECR. We do this when new images are uploaded with the same tag to save space.

To add the new image version
1. Go back to the ECR, and find the URI for the newest version of the image that you want to work on.
2. From the SageMaker Studio Control Panel, find the attached image that you'd like to revise at the bottom of the window.
3. Click on "Attach version" next to the image.

    <img src="https://user-images.githubusercontent.com/8322751/112722273-e9b88500-8ede-11eb-9279-988c5a9e02d0.png" width="600px">

4. Select "New Image Version"
5. Copy the URI into the textbox

    <img src="https://user-images.githubusercontent.com/8322751/112722358-69465400-8edf-11eb-9f08-8bf1af4de2a7.png" width="600px">

6. Click "Next"
7. Don't change anything on the Image properties or Studio configuration pages (unless the new image has new configurations you need to update)
8. Submit the new version
9. You should now see the "Latest version attached" number increase by 1 on the "Custom images attached to domain" box.
