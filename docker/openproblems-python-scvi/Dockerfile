FROM singlecellopenproblems/openproblems:latest

ARG NB_USER="sagemaker-user"
ARG NB_UID="1000"
ARG NB_GID="100"

USER root
WORKDIR /

# Install Python packages
COPY ./docker/openproblems-python-scvi/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

USER $NB_UID
WORKDIR /home/$NB_USER
