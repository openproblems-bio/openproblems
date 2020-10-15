FROM singlecellopenproblems/openproblems-r-base:latest

# Install R packages
COPY ./docker/openproblems-r-extras/r_requirements.txt ./r_requirements.txt
RUN for pkg in $(cat r_requirements.txt); do R -e "BiocManager::install('${pkg}')"; done

# Install Python packages
COPY ./docker/openproblems-r-extras/requirements.txt ./requirements.txt
RUN pip install -y -r requirements.txt

RUN apt-get clean -y && apt-get autoremove -y
