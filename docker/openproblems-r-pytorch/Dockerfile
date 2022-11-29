FROM singlecellopenproblems/openproblems-r-base:latest


ARG NB_USER="sagemaker-user"
ARG NB_UID="1000"
ARG NB_GID="100"

USER root
WORKDIR /

# Install R packages
COPY ./docker/openproblems-r-pytorch/r_requirements.txt ./r_requirements.txt
RUN R -e "source(\"install_renv.R\"); install_renv(\"r_requirements.txt\")"

# Install Python packages
COPY ./docker/openproblems-r-pytorch/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt


USER $NB_UID
WORKDIR /home/$NB_USER
