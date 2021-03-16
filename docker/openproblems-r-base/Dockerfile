FROM singlecellopenproblems/openproblems:latest

ARG NB_USER="sagemaker-user"
ARG NB_UID="1000"
ARG NB_GID="100"

USER $NB_UID

# Make the default shell bash (vs "sh") for a better Jupyter terminal UX
ENV SHELL=/bin/bash \
    NB_USER=$NB_USER \
    NB_UID=$NB_UID \
    NB_GID=$NB_GID \
    HOME=/home/$NB_USER \
    MINICONDA_VERSION=4.6.14 \
    CONDA_VERSION=4.6.14 \
    MINICONDA_MD5=718259965f234088d785cad1fbd7de03 \
    CONDA_DIR=/opt/conda \
    PATH=$CONDA_DIR/bin:${PATH}


USER root

# Install R
# RUN apt-get update -qq
# RUN apt-get -yq --no-install-suggests --no-install-recommends --allow-unauthenticated install dirmngr ca-certificates gnupg gpgv gfortran libblas-dev liblapack-dev
# RUN gpg --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
# RUN gpg --output key.gpg --export 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
# RUN apt-key add key.gpg
# RUN echo 'deb http://cloud.r-project.org/bin/linux/debian buster-cran40/' >> /etc/apt/sources.list
# RUN apt-get update -qq
# RUN apt-get install -yq --no-install-suggests --no-install-recommends r-base-dev=4.0\*
# RUN apt-get clean -y && apt-get autoremove -y

# Heavily inspired from https://github.com/jupyter/docker-stacks/blob/master/r-notebook/Dockerfile
# R system library pre-requisites
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    fonts-dejavu \
    unixodbc \
    unixodbc-dev \
    r-cran-rodbc \
    gfortran \
    gcc && \
    rm -rf /var/lib/apt/lists/* && \
    mkdir -p $CONDA_DIR && \
    chown -R $NB_USER:$NB_GID $CONDA_DIR && \
    # Fix for devtools https://github.com/conda-forge/r-devtools-feedstock/issues/4
    ln -s /bin/tar /bin/gtar

USER $NB_UID

# Install R packages
#RUN R -e "install.packages('renv'); renv::consent(TRUE)"
#COPY ./docker/openproblems-r-base/r_requirements.txt ./r_requirements.txt
# Install r_requirements.txt
#RUN for pkg in $(cat r_requirements.txt); do R -e "renv::install('${pkg}')"; done


ENV PATH=$CONDA_DIR/bin:${PATH}

# Install conda via Miniconda
RUN cd /tmp && \
    curl --silent --show-error --output miniconda-installer.sh https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh && \
    echo "${MINICONDA_MD5} *miniconda-installer.sh" | md5sum -c - && \
    /bin/bash miniconda-installer.sh -f -b -p $CONDA_DIR && \
    rm miniconda-installer.sh && \
    conda config --system --prepend channels conda-forge && \
    conda config --system --set auto_update_conda false && \
    conda config --system --set show_channel_urls true && \
    conda install --quiet --yes conda="${CONDA_VERSION%.*}.*" && \
    conda update --all --quiet --yes && \
    conda clean --all -f -y && \
    rm -rf /home/$NB_USER/.cache/yarn


# R packages and Python packages that are usable via "reticulate".
RUN conda install --quiet --yes \
    'r-base=4.0.0' \
    'r-caret=6.*' \
    'r-crayon=1.3*' \
    'r-devtools=2.3*' \
    'r-forecast=8.12*' \
    'r-hexbin=1.28*' \
    'r-htmltools=0.4*' \
    'r-htmlwidgets=1.5*' \
    'r-irkernel=1.1*' \
    'r-rmarkdown=2.2*' \
    'r-rodbc=1.3*' \
    'r-rsqlite=2.2*' \
    'r-shiny=1.4*' \
    'r-tidyverse=1.3*' \
    'unixodbc=2.3.*' \
    'r-tidymodels=0.1*' \
    'r-reticulate=1.*' \
    'RcppAnnoy'
    'bioc::scran'
    'IRkernel'
    && \
    pip install --quiet --no-cache-dir \
    'rpy2' \
    'anndata2ri>=1.0.6' \
    'boto3>1.0<2.0' \
    'sagemaker>2.0<3.0' && \
    conda clean --all -f -y

WORKDIR $HOME

USER $NB_UID
