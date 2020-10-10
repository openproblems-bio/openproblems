import os
import tempfile

import numpy as np
import anndata
import scprep
import scanpy as sc

from .utils import loader


@loader
def load_10xbrain(test=False):
    """Downloads 10x brain 1.3 million cell data from figshare"""
    URL = "https://ndownloader.figshare.com/files/25020254"

    with tempfile.TemporaryDirectory() as tempdir:
        filepath = os.path.join(tempdir, "10xbrain.h5ad")
        scprep.io.download.download_url(URL, filepath)
        adata = anndata.read_h5ad(filepath)

    if test:
        sc.pp.subsample(adata, n_obs=100)
        adata = adata[:, :100]
    return adata
