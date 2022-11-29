FROM singlecellopenproblems/openproblems:latest

ARG NB_USER="sagemaker-user"
ARG NB_UID="1000"
ARG NB_GID="100"

USER root
WORKDIR /

# Install pybedtools dependency
ARG BUILD_PACKAGES=""
ARG PACKAGE_VERSION=2.27.1
RUN apt-get update && \
    apt-get install --yes git openssl build-essential zlib1g-dev && \
    cd /tmp && \
    git clone https://github.com/arq5x/bedtools2.git && \
    cd bedtools2 && \
    git checkout v$PACKAGE_VERSION && \
    make && \
    mv bin/* /usr/local/bin && \
    cd /

# install dependencies and openproblems
COPY ./docker/openproblems-python-bedtools/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

USER $NB_UID
WORKDIR /home/$NB_USER
