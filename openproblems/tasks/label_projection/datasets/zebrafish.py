import os
import tempfile

import numpy as np
import anndata
import scprep
import scanpy as sc

from .utils import loader

URL = (
    "https://ndownloader.figshare.com/files/23992451?private_link=e3921450ec1bd0587870"
)


@loader
def load_zebrafish(test=False):
    if test:
        # load full data first, cached if available
        adata = load_zebrafish(test=False)
        sc.pp.subsample(adata, n_obs=100)
        adata = adata[:, :100]
        return adata
    else:
        with tempfile.TemporaryDirectory() as tempdir:
            filepath = os.path.join(tempdir, "zebrafish.h5ad")
            scprep.io.download.download_url(URL, filepath)
            adata = anndata.read_h5ad(filepath)
        return adata
