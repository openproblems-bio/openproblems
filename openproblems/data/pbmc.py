from . import utils

import os
import scanpy as sc
import scprep
import tempfile

URL = "https://ndownloader.figshare.com/files/24974582"


@utils.loader(data_url=URL)
def load_pbmc(test=False):
    """Download PBMC data from Figshare."""
    with tempfile.TemporaryDirectory() as tempdir:
        filepath = os.path.join(tempdir, "pbmc.h5ad")
        scprep.io.download.download_url(URL, filepath)
        adata = sc.read_h5ad(filepath)
        utils.filter_genes_cells(adata)

    if test:
        sc.pp.subsample(adata, n_obs=100)
        adata = adata[:, :1000]
        utils.filter_genes_cells(adata)

    return adata
