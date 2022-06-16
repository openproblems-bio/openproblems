FROM singlecellopenproblems/openproblems-r-base:latest

USER root
WORKDIR /

ARG NB_USER="sagemaker-user"
ARG NB_UID="1000"
ARG NB_GID="100"

RUN sed -i '$ d' /etc/apt/sources.list
RUN \
apt-get update --allow-releaseinfo-change && \
apt-get -y install --no-install-recommends gcc git python3-llvmlite && \
apt-get autoremove -y && \
rm -rf /var/lib/apt/lists/*

# Install Python packages
COPY ./docker/openproblems-python-batch-integration/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt --no-binary annoy

USER $NB_UID
WORKDIR /home/$NB_USER
