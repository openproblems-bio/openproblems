FROM singlecellopenproblems/openproblems-r-base:latest


ARG NB_USER="sagemaker-user"
ARG NB_UID="1000"
ARG NB_GID="100"

USER root
WORKDIR /

RUN apt-get update && apt-get install -y \
    libhdf5-dev hdf5-tools libgeos-dev \
    libharfbuzz-dev libfribidi-dev

RUN apt-get clean autoclean && \
    apt-get autoremove --yes && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/ && \
    apt-get install -y libhdf5-dev=1.10.6+repack-4+deb11u1

# Install R packages
COPY ./docker/openproblems-r-extras/r_requirements.txt ./r_requirements.txt
RUN R -e "source(\"install_renv.R\"); install_renv(\"r_requirements.txt\")"

# install dependencies and openproblems
COPY ./docker/openproblems-r-extras/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt --editable /usr/src/singlecellopenproblems[r]

# Fix permissions
RUN chown -R $NB_USER:$NB_GID /home/$NB_USER

USER $NB_UID
WORKDIR /home/$NB_USER
