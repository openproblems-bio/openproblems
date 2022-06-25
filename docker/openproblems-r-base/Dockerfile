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
ENV R_VERSION="4.2.1"
RUN apt-get update -qq
RUN apt-get -yq --no-install-suggests --no-install-recommends --allow-unauthenticated install \
    dirmngr ca-certificates gnupg gpgv gfortran liblapack-dev gdebi-core \
    libopenblas-dev libopenblas-pthread-dev libopenblas0 libopenblas0-pthread libpaper-utils libpaper1 zip
RUN apt-get clean -y && apt-get autoremove -y
RUN curl https://cdn.rstudio.com/r/debian-11/pkgs/r-${R_VERSION}_1_amd64.deb --output r-${R_VERSION}_3_amd64.deb
RUN gdebi -n r-${R_VERSION}_3_amd64.deb
RUN update-alternatives --install /usr/bin/R R /opt/R/${R_VERSION}/bin/R 1
RUN update-alternatives --install /usr/bin/Rscript Rscript /opt/R/${R_VERSION}/bin/Rscript 1
ENV R_HOME=/opt/R/${R_VERSION}/lib/R
RUN mkdir -p ${R_HOME}/etc
RUN bash -c "echo \"local({r <- getOption('repos'); r['CRAN'] <- 'https://cloud.r-project.org/'; options(repos = r)})\" >> ${R_HOME}/etc/Rprofile.site"

# Install R packages
RUN R -e "install.packages('renv'); renv::consent(TRUE)"
COPY ./docker/openproblems-r-base/r_requirements.txt ./r_requirements.txt
RUN R -e "renv::install(scan('r_requirements.txt', sep='\n', what=character()))"

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
