FROM singlecellopenproblems/openproblems-r-base:latest


ARG NB_USER="sagemaker-user"
ARG NB_UID="1000"
ARG NB_GID="100"

USER root
WORKDIR /

RUN apt-get update && apt-get install -y libhdf5-dev hdf5-tools \
&& apt-get clean autoclean && apt-get autoremove --yes && rm -rf /var/lib/{apt,dpkg,cache,log}/
RUN apt-get install -y libhdf5-dev=1.10.6+repack-4+deb11u1

# Install R packages
COPY ./docker/openproblems-r-extras/r_requirements.txt ./r_requirements.txt
RUN \
    R -e "renv::install(scan('r_requirements.txt', sep='\n', what=character()))" \
    && R -e 'devtools::install_github("immunogenomics/harmony")' \
    && R -e 'devtools::install_github("welch-lab/liger@v0.5.0")'



# Install Python packages
COPY ./docker/openproblems-r-extras/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt


USER $NB_UID
WORKDIR /home/$NB_USER
