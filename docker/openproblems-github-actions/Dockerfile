FROM singlecellopenproblems/openproblems-r-base:latest

USER root
WORKDIR /

# Install Docker
RUN apt-get update
RUN apt-get install -y ca-certificates curl gnupg lsb-release
RUN sh -c 'curl -fsSL https://download.docker.com/linux/ubuntu/gpg | gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg'
RUN sh -c 'echo \
	"deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/debian \
	$(lsb_release -cs) stable" | tee /etc/apt/sources.list.d/docker.list'
RUN apt-get update
RUN apt-get install -y docker-ce docker-ce-cli containerd.io

# install Python packages
COPY ./docker/openproblems-github-actions/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Install R packages
COPY ./docker/openproblems-github-actions/r_requirements.txt ./r_requirements.txt
RUN R -e "source(\"install_renv.R\"); install_renv(\"r_requirements.txt\")"

# Install Python packages
RUN pip install --no-cache-dir -U /usr/src/singlecellopenproblems[test,r,evaluate]
