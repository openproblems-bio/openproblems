import scanpy as sc
import scprep
import tempfile
import os
from . import utils

URL = "https://ndownloader.figshare.com/files/24130931"


@utils.loader
def load_senis_muris(test=False):
    """Download Tabula Senis."""
    if test:
        adata = load_senis_muris(test=False)
        sc.pp.subsample(adata, fraction=0.1)
        utils.filter_genes_cells(adata)
        return adata
    else:
        with tempfile.TemporaryDirectory() as tempdir:
            filepath = os.path.join(tempdir, "senis_muris.h5ad")
            scprep.io.download.download_url(URL, filepath)
            adata = sc.read(filepath)
            utils.filter_genes_cells(adata)
            return adata
