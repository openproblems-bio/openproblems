FROM singlecellopenproblems/openproblems-r-base:latest

# Install R packages
COPY ./docker/openproblems-r-extras/r_requirements.txt ./r_requirements.txt
RUN renv::install(scan('r_requirements.txt', sep='\n', what=character()))

# Install Python packages
COPY ./docker/openproblems-r-extras/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt
