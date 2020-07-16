import os
import tempfile

import numpy as np
import anndata
import scprep
import scanpy as sc

from .utils import loader

URL = (
    "https://ndownloader.figshare.com/files/23820602?private_link=a23de5c2117a69d64af4"
)


@loader
def load_zebrafish(test=False):

    with tempfile.TemporaryDirectory() as tempdir:
        filepath = os.path.join(tempdir, "zebrafish.h5ad")
        scprep.io.download.download_url(URL, filepath)
        adata = anndata.read_h5ad(filepath)

    if test:
        sc.pp.subsample(adata, n_obs=100)
    return adata
