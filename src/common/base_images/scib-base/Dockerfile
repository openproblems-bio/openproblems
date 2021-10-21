# version 0.1.2
FROM mumichae/scanpy-r-micromamba:0.1.2
COPY env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml
RUN git clone https://github.com/theislab/scib.git
RUN cd scib && \
    pip install .
RUN micromamba clean --all --yes
