FROM rackspacedot/python37:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN pip install --no-cache-dir -U pip wheel setuptools cmake

# Install R
RUN apt install dirmngr apt-transport-https ca-certificates software-properties-common gnupg2
RUN apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/debian stretch-cran35/'
RUN apt-get update -qq
RUN apt-get install r-base-core=3.6\* -y
RUN export R_LIBS_USER="$HOME/R/Library"
RUN echo ".libPaths(c('$R_LIBS_USER', .libPaths()))" >> $HOME/.Rprofile

# Install single-cell open problems
RUN pip install git+https://github.com/singlecellopenproblems/SingleCellOpenProblems.git#egg=openproblems[methods]

RUN apt-get clean -y && apt-get autoremove -y
