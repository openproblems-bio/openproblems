FROM singlecellopenproblems/openproblems-r-base:latest

ARG NB_USER="sagemaker-user"
ARG NB_UID="1000"
ARG NB_GID="100"

USER root
WORKDIR /

# Install R packages
COPY ./docker/openproblems-r-extras/r_requirements.txt ./r_requirements.txt
RUN sudo R -e "renv::install(scan('r_requirements.txt', sep='\n', what=character()))"

# Install Python packages
COPY ./docker/openproblems-r-extras/requirements.txt ./requirements.txt
RUN sudo pip install --no-cache-dir -r requirements.txt

USER $NB_UID
WORKDIR /home/$NB_USER
