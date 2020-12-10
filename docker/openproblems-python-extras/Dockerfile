FROM singlecellopenproblems/openproblems:latest

# Install pybedtools dependency
apt install bedtools

# Install Python packages
COPY ./docker/openproblems-python-extras/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt
