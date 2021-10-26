# version 1.0.0
FROM mumichae/scanpy-r-micromamba:0.1.2
COPY env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml
RUN micromamba clean --all --yes
