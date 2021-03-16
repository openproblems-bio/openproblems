FROM python:3.8-buster

# Removing on 3/16 to see if essential
# ENV DEBIAN_FRONTEND=noninteractive

# Setting up Sagemaker Studio Image from example
# https://github.com/aws-samples/sagemaker-studio-custom-image-samples/blob/main/examples/echo-kernel-image/Dockerfile

ARG NB_USER="sagemaker-user"
ARG NB_UID="1000"
ARG NB_GID="100"

# Setup the "sagemaker-user" user with root privileges.
RUN \
   apt-get update && \
   apt-get install -y sudo && \
   useradd -m -s /bin/bash -N -u $NB_UID $NB_USER && \
   chmod g+w /etc/passwd && \
   echo "${NB_USER}    ALL=(ALL)    NOPASSWD:    ALL" >> /etc/sudoers && \
   # Prevent apt-get cache from being persisted to this layer.
   rm -rf /var/lib/apt/lists/*

# update pip
RUN python3 -m pip install --no-cache-dir -U pip

# install AWS CLI
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install

# install dependencies
COPY ./docker/openproblems/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Install single-cell open problems
COPY . /usr/src/singlecellopenproblems
RUN cd /usr/src/singlecellopenproblems && git clean -fxd
RUN pip install --no-cache-dir --editable /usr/src/singlecellopenproblems

# Make the default shell bash (vs "sh") for a better Jupyter terminal UX
ENV SHELL=/bin/bash

USER $NB_UID
