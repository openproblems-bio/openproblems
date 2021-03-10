FROM singlecellopenproblems/openproblems-r-base:latest

RUN apt-get update && apt-get -y install gcc

# Install Python packages
COPY ./docker/openproblems-python-batch-integration/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt
