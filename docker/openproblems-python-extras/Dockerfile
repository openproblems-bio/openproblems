FROM singlecellopenproblems/openproblems:latest

# Install Python packages
COPY ./docker/openproblems-python-extras/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -U -r requirements.txt
