FROM singlecellopenproblems/openproblems:latest

# Install R
RUN apt-get update -qq
RUN apt-get -yq --no-install-suggests --no-install-recommends --allow-unauthenticated install dirmngr ca-certificates gnupg gpgv gfortran libblas-dev liblapack-dev
RUN gpg --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
RUN gpg --output key.gpg --export 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
RUN apt-key add key.gpg
RUN echo 'deb http://cloud.r-project.org/bin/linux/debian buster-cran40/' >> /etc/apt/sources.list
RUN apt-get update -qq
RUN apt-get install -yq --no-install-suggests --no-install-recommends r-base-dev=4.0\*
RUN apt-get clean -y && apt-get autoremove -y

# Install R packages
RUN R -e "install.packages('renv'); renv::consent(TRUE)"
COPY ./docker/openproblems-r-base/r_requirements.txt ./r_requirements.txt
# Install r_requirements.txt
RUN for pkg in $(cat r_requirements.txt); do R -e "renv::install('${pkg}')"; done

# Install Python packages
COPY ./docker/openproblems-r-base/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -U -r requirements.txt
