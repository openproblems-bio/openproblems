FROM singlecellopenproblems/openproblems-r-base:latest

RUN apt-get update
RUN apt-get install -y libhdf5-dev

# Install R packages
COPY ./docker/openproblems-r-extras/r_requirements.txt ./r_requirements.txt
RUN R -e "renv::install(scan('r_requirements.txt', sep='\n', what=character()))"
RUN R -e 'install.packages("devtools")'
RUN R -e 'devtools::install_github("immunogenomics/harmony")'
RUN R -e 'devtools::install_github("welch-lab/liger@v0.5.0")'


# Install Python packages
COPY ./docker/openproblems-r-extras/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# RUN R -e 'BiocManager::install("SingleCellExperiment")'
