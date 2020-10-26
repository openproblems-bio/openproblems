FROM singlecellopenproblems/openproblems-r-base:latest

# Install R packages
COPY ./docker/openproblems-r-extras/r_requirements.txt ./r_requirements.txt
RUN for pkg in $(cat r_requirements.txt); do R -e "BiocManager::install('${pkg}')" 2>&1 > /dev/null | grep -v "^\*\*"; done

# Install Python packages
COPY ./docker/openproblems-r-extras/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -U -r requirements.txt
