import os
import tempfile

import numpy as np
import anndata
import scprep
import scanpy as sc

from .utils import loader


@loader
def load_zebrafish(test=False):
    """Downloads zebrafish data from figshare"""
    URL = "https://ndownloader.figshare.com/files/23992451?private_link=e3921450ec1bd0587870"

    with tempfile.TemporaryDirectory() as tempdir:
        filepath = os.path.join(tempdir, "zebrafish.h5ad")
        scprep.io.download.download_url(URL, filepath)
        adata = anndata.read_h5ad(filepath)

    if test:
        sc.pp.subsample(adata, n_obs=100)
        adata = adata[:, :100]
    return adata
