FROM rackspacedot/python37:latest

ENV DEBIAN_FRONTEND=noninteractive

COPY ./docker/openproblems/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Install single-cell open problems
COPY . /usr/src/singlecellopenproblems
RUN pip install --no-cache-dir /usr/src/singlecellopenproblems
