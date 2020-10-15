FROM rackspacedot/python37:latest

ENV DEBIAN_FRONTEND=noninteractive

ARG SOURCEDIR=.

RUN pip install --no-cache-dir -U pip wheel setuptools cmake

# Install single-cell open problems
COPY $SOURCEDIR /usr/src/singlecellopenproblems
RUN pip install -e /usr/src/singlecellopenproblems

RUN apt-get clean -y && apt-get autoremove -y
