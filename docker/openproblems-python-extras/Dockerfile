FROM singlecellopenproblems/openproblems:latest

# Install pybedtools dependency
ARG BUILD_PACKAGES=""
ARG PACKAGE_VERSION=2.27.1
RUN apt-get update && \
    apt-get install --yes git openssl build-essential zlib1g-dev && \
    cd /tmp && \
    git clone https://github.com/arq5x/bedtools2.git && \
    cd bedtools2 && \
    git checkout v$PACKAGE_VERSION && \
    make && \
    mv bin/* /usr/local/bin && \
    cd /

# Install Python packages
COPY ./docker/openproblems-python-extras/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt
