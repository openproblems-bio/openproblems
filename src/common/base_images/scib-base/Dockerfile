# version 1.0.2
FROM mumichae/scanpy-r-micromamba:1.9.1
COPY env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml
RUN Rscript -e "devtools::install_github('theislab/kBET')"
RUN micromamba clean --all --yes
