FROM debian:10

ENV DEBIAN_FRONTEND=noninteractive

# Install system libraries required for python and R installations
RUN apt-get update && apt-get install -y --no-install-recommends build-essential apt-utils ca-certificates zlib1g-dev gfortran locales libxml2-dev libcurl4-openssl-dev libssl-dev libzmq3-dev libreadline6-dev xorg-dev libcairo-dev libpango1.0-dev libbz2-dev liblzma-dev libffi-dev libsqlite3-dev nodejs npm

# Install common linux tools
RUN apt-get update && apt-get install -y --no-install-recommends wget curl htop less nano vim emacs git

# Configure default locale
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8

# Download and compile python from source
WORKDIR /opt/python
RUN wget https://www.python.org/ftp/python/3.7.7/Python-3.7.7.tgz \
    && tar zxfv Python-3.7.7.tgz && rm Python-3.7.7.tgz

WORKDIR /opt/python/Python-3.7.7
RUN ./configure --enable-optimizations --with-lto --prefix=/opt/python/ \
    && make \
    && make install

WORKDIR /opt/python
RUN rm -rf /opt/python/Python-3.7.7 \
    && ln -s /opt/python/bin/python3 /opt/python/bin/python \
    && ln -s /opt/python/bin/pip3 /opt/python/bin/pip
ENV PATH="/opt/python/bin:${PATH}"

RUN pip install --no-cache-dir -U pip wheel setuptools cmake

# Install single-cell open problems
RUN pip install git+https://github.com/singlecellopenproblems/SingleCellOpenProblems.git

RUN apt-get clean -y && apt-get autoremove -y

