FROM openproblems-python-data-integration:latest

# Install Python packages
COPY ./. ./.
copy ./docker/openproblems-r-custom/requirements.txt ./requirements.txt
RUN pip install --no-cache-dir -e .
RUN pip install --no-cache-dir -r ./requirements.txt
