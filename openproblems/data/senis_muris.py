import os
import tempfile


import scprep
import scanpy as sc

from .utils import loader


@loader
def load_senis_muris(test=False):
    url = "https://ndownloader.figshare.com/files/24130931"
    with tempfile.TemporaryDirectory() as tempdir:
        filepath = os.path.join(tempdir, "senis_muris_sub.h5ad")
        scprep.io.download.download_url(url, filepath)
        adata = sc.read_h5ad(filepath)
    if test:
        sc.pp.subsample(adata, fraction=0.1)
    return adata
