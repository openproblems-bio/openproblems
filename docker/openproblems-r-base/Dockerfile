FROM singlecellopenproblems/openproblems:latest

# Install R
RUN apt install dirmngr ca-certificates gnupg gpgv
RUN apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
RUN echo 'deb http://cloud.r-project.org/bin/linux/debian stretch-cran35/' >> /etc/apt/sources.list
RUN apt-get update -qq
RUN apt-get install r-base-core=3.6\* -y

# Install R packages
RUN R -e "install.packages('BiocManager'); BiocManager::install('scran')"
RUN R -e "BiocManager::install('batchelor')"

# Install Python packages
RUN pip install git+https://github.com/singlecellopenproblems/SingleCellOpenProblems.git#egg=openproblems[r]


RUN apt-get clean -y && apt-get autoremove -y
