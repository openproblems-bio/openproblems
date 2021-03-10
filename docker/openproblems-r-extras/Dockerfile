FROM singlecellopenproblems/openproblems-r-base:latest

RUN apt-get update
RUN apt-get install -y libhdf5-dev

# Install R packages
COPY ./docker/openproblems-r-extras/r_requirements.txt ./r_requirements.txt
RUN for pkg in $(cat r_requirements.txt); do R -e "renv::install('${pkg}')"; done
RUN R -e 'install.packages("devtools")'
RUN R -e 'devtools::install_github("immunogenomics/harmony")'

# Install Python packages
COPY ./docker/openproblems-r-extras/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

RUN R -e 'BiocManager::install("SingleCellExperiment")'
