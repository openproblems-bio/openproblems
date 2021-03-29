FROM singlecellopenproblems/openproblems-r-base:latest


ARG NB_USER="sagemaker-user"
ARG NB_UID="1000"
ARG NB_GID="100"

USER root
WORKDIR /

RUN sudo apt-get update && sudo apt-get install -y libhdf5-dev hdf5-tools

# Install R packages
COPY ./docker/openproblems-r-extras/r_requirements.txt ./r_requirements.txt
RUN sudo R -e "renv::install(scan('r_requirements.txt', sep='\n', what=character()))"
RUN sudo R -e 'install.packages("devtools")'
RUN sudo R -e 'devtools::install_github("immunogenomics/harmony")'
RUN sudo R -e 'devtools::install_github("welch-lab/liger@v0.5.0")'


# Install Python packages
COPY ./docker/openproblems-r-extras/requirements.txt ./requirements.txt
RUN sudo pip install --no-cache-dir -r requirements.txt


USER $NB_UID
WORKDIR /home/$NB_USER
