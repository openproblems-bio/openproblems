FROM rackspacedot/python38:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN python3 -m pip install --no-cache-dir -U pip

COPY ./docker/openproblems/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Install single-cell open problems
COPY . /usr/src/singlecellopenproblems
RUN cd /usr/src/singlecellopenproblems && git clean -fxd
RUN pip install --no-cache-dir /usr/src/singlecellopenproblems
