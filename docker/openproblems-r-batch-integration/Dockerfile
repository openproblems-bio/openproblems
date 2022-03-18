FROM singlecellopenproblems/openproblems-python37-scgen:latest

ENV DEBIAN_FRONTEND=noninteractive
ARG NB_USER="sagemaker-user"
ARG NB_UID="1000"
ARG NB_GID="100"
USER root



RUN mamba install -y -c bioconda -c conda-forge bioconductor-batchelor bioconductor-complexheatmap r-seurat r-dplyr r-devtools r-harmony 
COPY ./docker/openproblems-r-batch-integration/requirements_mamba.txt requirements_mamba.txt
RUN mamba install -y --file requirements_mamba.txt
RUN R -e 'devtools::install_github("welch-lab/liger@v0.5.0")'

# Install openproblems
COPY . /usr/src/singlecellopenproblems
RUN cd /usr/src/singlecellopenproblems && git clean -fxdq
RUN pip install --no-cache-dir --editable /usr/src/singlecellopenproblems

# Overwrite kernel.json to use system Python install
COPY ./docker/openproblems/kernelspec.json /usr/local/share/jupyter/kernels/python3/kernel.json

# Make the default shell bash (vs "sh") for a better Jupyter terminal UX
ENV SHELL=/bin/bash

USER $NB_UID
WORKDIR /home/$NB_USER

