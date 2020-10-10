import os
import tempfile

import numpy as np
import scprep
import scanpy as sc

from .utils import loader
from ..tools import normalize as n


@loader
def load_pbmc(test=False):
    """Downloads PBMC data from figshare"""
    URL = "https://ndownloader.figshare.com/files/24974582"

    with tempfile.TemporaryDirectory() as tempdir:
        filepath = os.path.join(tempdir, "pbmc.h5ad")
        scprep.io.download.download_url(URL, filepath)
        adata = sc.read_h5ad(filepath)

    if test:
        sc.pp.subsample(adata, n_obs=100)
        adata = adata[:, :1000]

    return adata
