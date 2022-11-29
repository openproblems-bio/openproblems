FROM singlecellopenproblems/openproblems:latest

# Install Python packages
COPY ./docker/openproblems-python-scvi/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -r requirements.txt
