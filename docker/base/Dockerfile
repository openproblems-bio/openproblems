FROM rackspacedot/python37:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN pip install --no-cache-dir -U pip wheel setuptools cmake

# Install R
RUN apt install dirmngr ca-certificates gnupg gpgv
RUN apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
RUN echo 'deb http://cloud.r-project.org/bin/linux/debian stretch-cran35/' >> /etc/apt/sources.list
RUN apt-get update -qq
RUN apt-get install r-base-core=3.6\* -y

# R local library setup
RUN echo 'export R_LIBS_USER="$HOME/R/Library"' >> $HOME/.bashrc
RUN export R_LIBS_USER="$HOME/R/Library"
RUN echo ".libPaths(c('$R_LIBS_USER', .libPaths()))" >> $HOME/.Rprofile

# Install single-cell open problems
RUN pip install git+https://github.com/singlecellopenproblems/SingleCellOpenProblems.git#egg=openproblems[methods]

RUN apt-get clean -y && apt-get autoremove -y
