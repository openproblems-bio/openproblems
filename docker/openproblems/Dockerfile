FROM rackspacedot/python38:latest

ENV DEBIAN_FRONTEND=noninteractive

# update pip
RUN python3 -m pip install --no-cache-dir -U pip

# install AWS CLI
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip -q awscliv2.zip
RUN ./aws/install

# install dependencies
COPY ./docker/openproblems/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Install single-cell open problems
COPY . /usr/src/singlecellopenproblems
RUN cd /usr/src/singlecellopenproblems && git clean -fxd
RUN pip install --no-cache-dir --editable /usr/src/singlecellopenproblems
