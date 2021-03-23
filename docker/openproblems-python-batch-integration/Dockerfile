FROM singlecellopenproblems/openproblems-r-base:latest

RUN apt-get update && apt-get -y install gcc git python3-llvmlite

# Install Python packages
COPY ./docker/openproblems-python-batch-integration/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt
RUN pip install --no-cache-dir --no-deps louvain==0.6.1 umap-learn==0.4.6 rpy2==3.4.2


