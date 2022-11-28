FROM singlecellopenproblems/openproblems-r-base:latest

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

# Install Python packages
COPY ./docker/openproblems-r-extras/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt
