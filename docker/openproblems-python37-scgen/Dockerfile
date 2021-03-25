FROM python:3.7

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
RUN apt-get -yq --no-install-suggests --no-install-recommends --allow-unauthenticated install dirmngr ca-certificates gnupg gpgv gfortran libblas-dev liblapack-dev
RUN gpg --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
RUN gpg --output key.gpg --export 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
RUN apt-key add key.gpg
RUN echo 'deb http://cloud.r-project.org/bin/linux/debian buster-cran40/' >> /etc/apt/sources.list
RUN apt-get update -qq
RUN apt-get install -yq --no-install-suggests --no-install-recommends r-base-dev=4.0\*
RUN apt-get clean -y && apt-get autoremove -y

# update pip
RUN python3 -m pip install --no-cache-dir -U pip

# install AWS CLI
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip -q awscliv2.zip
RUN ./aws/install

# install dependencies
COPY ./docker/openproblems/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Install Saucie
RUN git clone https://github.com/KrishnaswamyLab/SAUCIE
RUN mv SAUCIE /usr/local/lib/python3.7/site-packages/

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
