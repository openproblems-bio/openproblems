FROM singlecellopenproblems/openproblems-r-base:latest

RUN sed -i '$ d' /etc/apt/sources.list
RUN \
apt-get update --allow-releaseinfo-change && \
apt-get -y install --no-install-recommends gcc git python3-llvmlite && \
apt-get autoremove -y && \
rm -rf /var/lib/apt/lists/*

# Install Python packages
COPY ./docker/openproblems-python-batch-integration/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt
# force reinstall annoy addresses https://github.com/spotify/annoy/issues/513
RUN pip install --no-cache-dir --force annoy==1.17.0
