FROM singlecellopenproblems/openproblems-r-base:latest

# Install R packages
COPY ./docker/openproblems-r-extras/r_requirements.txt ./r_requirements.txt
RUN for pkg in $(cat r_requirements.txt); do R -e "renv::install('${pkg}')"; done

# Install Python packages
COPY ./docker/openproblems-r-extras/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt
