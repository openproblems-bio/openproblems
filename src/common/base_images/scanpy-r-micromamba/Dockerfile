FROM mambaorg/micromamba:0.14.0
COPY env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes
RUN pip install anndata2ri==1.0.6
WORKDIR /home/micromamba
