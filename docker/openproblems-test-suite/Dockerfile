FROM singlecellopenproblems/openproblems-r-base:latest

# Need system python 3.5
RUN update-alternatives --remove python /usr/local/bin/python
RUN update-alternatives --install /usr/local/bin/python python /usr/bin/python3.5 10

# Dependencies
RUN apt-get update
RUN apt-get install -y lsb-core software-properties-common apt-transport-https init

# Install Docker
RUN curl -fsSL https://download.docker.com/linux/debian/gpg | apt-key add -
RUN apt-key fingerprint 0EBFCD88
RUN add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/debian \
   $(lsb_release -cs) \
   stable"
RUN apt-get update
RUN apt-get install -y docker-ce docker-ce-cli containerd.io
RUN systemctl enable docker

# Install singularity
RUN wget -O- http://neuro.debian.net/lists/xenial.us-ca.full >> /etc/apt/sources.list.d/neurodebian.sources.list
RUN curl http://neuro.debian.net/_static/neuro.debian.net.asc | apt-key add
RUN apt-get update
RUN apt-get install -y singularity-container

RUN python3.7 -m pip install /usr/src/singlecellopenproblems[test]
RUN python3.7 -m pip uninstall -y openproblems
