FROM rackspacedot/python37:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN pip install --no-cache-dir -U pip wheel setuptools cmake

# Install R
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial-cran35/" | sudo tee -a /etc/apt/sources.list
sudo apt-get update -qq
sudo apt-get install r-base-core=3.6\* -y
export R_LIBS_USER="$HOME/R/Library"
echo ".libPaths(c('$R_LIBS_USER', .libPaths()))" >> $HOME/.Rprofile

# Install single-cell open problems
RUN pip install git+https://github.com/singlecellopenproblems/SingleCellOpenProblems.git#egg=openproblems[methods]

RUN apt-get clean -y && apt-get autoremove -y
