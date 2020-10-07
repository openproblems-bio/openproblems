FROM rackspacedot/python37:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN pip install --no-cache-dir -U pip wheel setuptools cmake

# Install R
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial-cran35/" >> /etc/apt/sources.list
RUN apt-get update -qq
RUN apt-get install r-base-core=3.6\* -y
RUN export R_LIBS_USER="$HOME/R/Library"
RUN echo ".libPaths(c('$R_LIBS_USER', .libPaths()))" >> $HOME/.Rprofile

# Install single-cell open problems
RUN pip install git+https://github.com/singlecellopenproblems/SingleCellOpenProblems.git#egg=openproblems[methods]

RUN apt-get clean -y && apt-get autoremove -y
