FROM singlecellopenproblems/openproblems-r-base:latest

RUN apt-get update && apt-get -y install gcc git python3-llvmlite

# Install Python packages
COPY ./docker/openproblems-python-batch-integration/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Install Saucie
RUN git clone https://github.com/KrishnaswamyLab/SAUCIE
RUN mv SAUCIE /usr/local/lib/python3.8/site-packages/
