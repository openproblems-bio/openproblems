FROM rackspacedot/python37:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN pip install --no-cache-dir -U pip wheel setuptools cmake

# Install single-cell open problems
RUN pip install git+https://github.com/singlecellopenproblems/SingleCellOpenProblems.git

RUN apt-get clean -y && apt-get autoremove -y
