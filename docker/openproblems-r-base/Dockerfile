FROM singlecellopenproblems/openproblems:latest

# Install R
RUN apt install dirmngr ca-certificates gnupg gpgv
RUN apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF'
RUN echo 'deb http://cloud.r-project.org/bin/linux/debian stretch-cran35/' >> /etc/apt/sources.list
RUN apt-get update -qq
RUN apt-get install r-base-core=3.6\* -y

# Install R packages
COPY ./docker/openproblems-r-base/r_requirements.txt ./r_requirements.txt
RUN R -e "install.packages('BiocManager')"
RUN for pkg in $(cat r_requirements.txt); do R -e "BiocManager::install('${pkg}')"; done

# Install Python packages
COPY . /usr/src/singlecellopenproblems
RUN pip install -e /usr/src/singlecellopenproblems[r]

RUN apt-get clean -y && apt-get autoremove -y
