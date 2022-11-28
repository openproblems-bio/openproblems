FROM singlecellopenproblems/openproblems-r-base:latest

# Install R packages
COPY ./docker/openproblems-r-pytorch/r_requirements.txt ./r_requirements.txt
RUN R -e "source(\"install_renv.R\"); install_renv(\"r_requirements.txt\")"

# Install Python packages
COPY ./docker/openproblems-r-pytorch/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt
