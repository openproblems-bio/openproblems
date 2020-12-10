FROM singlecellopenproblems/openproblems:latest

# Install pybedtools dependency
RUN apt-get install -yq --no-install-suggests --no-install-recommends bedtools


# Install Python packages
COPY ./docker/openproblems-python-extras/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt
