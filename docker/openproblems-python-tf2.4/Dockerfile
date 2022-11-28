FROM singlecellopenproblems/openproblems:latest

# Install Python packages
COPY ./docker/openproblems-python-tf2.4/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# tensorflow downgrades numpy and h5py (and therefore anndata)
RUN pip install --no-cache-dir -e /usr/src/singlecellopenproblems
