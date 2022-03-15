#FROM python:3.7
FROM condaforge/mambaforge

ENV DEBIAN_FRONTEND=noninteractive
ARG NB_USER="sagemaker-user"
ARG NB_UID="1000"
ARG NB_GID="100"

RUN \
    apt-get update && \
    apt-get install -qy sudo && \
    useradd -m -s /bin/bash -N -u $NB_UID $NB_USER && \
    chmod g+w /etc/passwd && \
    echo "${NB_USER}    ALL=(ALL)    NOPASSWD:    ALL" >> /etc/sudoers && \
    # Prevent apt-get cache from being persisted to this layer.
    find /var/lib/apt/lists/ ! -type d -exec rm '{}' \;

RUN apt-get update -qq
RUN apt-get clean -y && apt-get autoremove -y && apt-get install -y cmake
RUN apt-get install -y git
# RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'
# RUN echo 'deb http://cloud.r-project.org/bin/linux/debian bullseye-cran40/' >> /etc/apt/sources.list
# RUN apt-get update -qq
# RUN apt-get install -yq --no-install-suggests --no-install-recommends r-base-dev=4.1\*
# RUN apt-get clean -y && apt-get autoremove -y

RUN mamba install -y python=3.7
RUN mamba install -c bioconda -y bioconductor-scran

# update pip
RUN python3 -m pip install --no-cache-dir -U pip

# install AWS CLI
RUN apt-get install -y unzip
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip -q awscliv2.zip
RUN ./aws/install

# install dependencies
COPY ./docker/openproblems-python37-saucie/requirements.txt requirements.txt
RUN pip install --no-cache-dir git+https://github.com/theislab/scib
RUN pip install --no-cache-dir -r requirements.txt
RUN git clone https://github.com/KrishnaswamyLab/SAUCIE
RUN cp -r SAUCIE /opt/conda/lib/python3.7/site-packages/


# Install single-cell open problems
COPY . /usr/src/singlecellopenproblems
RUN cd /usr/src/singlecellopenproblems && git clean -fxdq
RUN pip install --no-cache-dir --editable /usr/src/singlecellopenproblems

# Overwrite kernel.json to use system Python install
COPY ./docker/openproblems/kernelspec.json /usr/local/share/jupyter/kernels/python3/kernel.json

# Make the default shell bash (vs "sh") for a better Jupyter terminal UX
ENV SHELL=/bin/bash

USER $NB_UID
WORKDIR /home/$NB_USER
