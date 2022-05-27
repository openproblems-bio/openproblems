# version 1.9.1
FROM mambaorg/micromamba:0.14.0
COPY env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes
WORKDIR /home/micromamba
