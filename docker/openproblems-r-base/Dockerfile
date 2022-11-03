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
    HOME=/home/$NB_USER

USER root
WORKDIR /

# Install R
RUN apt-get update -qq
RUN apt-get -yq --no-install-suggests --no-install-recommends --allow-unauthenticated install dirmngr ca-certificates gnupg gpgv gfortran libblas-dev liblapack-dev
RUN gpg --keyserver keyserver.ubuntu.com --recv-key '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'
RUN gpg --output cran41.gpg --export '95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'
RUN cp cran41.gpg /etc/apt/trusted.gpg.d/
RUN echo 'deb http://cloud.r-project.org/bin/linux/debian bullseye-cran40/' >> /etc/apt/sources.list
RUN apt-get update -qq
RUN apt-get install -yq --no-install-suggests --no-install-recommends r-base-dev=4.2\*
RUN apt-get clean -y && apt-get autoremove -y
ENV R_HOME=/usr/lib/R
ENV BIOCVERSION="3.15"

# Install R packages
RUN R -e "install.packages('renv'); renv::consent(TRUE)"
COPY ./scripts/install_renv.R ./
COPY ./docker/openproblems-r-base/r_requirements.txt ./r_requirements.txt
RUN R -e "source(\"install_renv.R\"); install_renv(\"r_requirements.txt\")"

# Install kernelspec
RUN R -e "IRkernel::installspec(user = FALSE);"

# Install single-cell open problems with R requirements
COPY . /usr/src/singlecellopenproblems
RUN cd /usr/src/singlecellopenproblems && git clean -fxdq
RUN pip install --no-cache-dir --editable /usr/src/singlecellopenproblems[r]

# Fix permissions
RUN chown -R $NB_USER:$NB_GID /home/$NB_USER

USER $NB_UID
WORKDIR $HOME
