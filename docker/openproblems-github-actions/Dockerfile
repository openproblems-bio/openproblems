FROM singlecellopenproblems/openproblems-r-base:latest

USER root
WORKDIR /

# Install Docker
RUN sudo apt-get update
RUN sudo apt-get install -y ca-certificates curl gnupg lsb-release
RUN sh -c 'curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg'
RUN sh -c 'echo \
	"deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/debian \
	$(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list'
RUN sudo apt-get update
RUN sudo apt-get install -y docker-ce docker-ce-cli containerd.io

# Install R packages
COPY ./.pre-commit.r_requirements.txt ./r_requirements.txt
RUN sudo R -e "renv::install(scan('r_requirements.txt', sep='\n', what=character()))"

# Install Python packages
RUN sudo pip install --no-cache-dir -U /usr/src/singlecellopenproblems[test,r,evaluate]
