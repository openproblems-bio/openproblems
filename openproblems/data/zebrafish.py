from . import utils

import anndata
import os
import scprep
import tempfile

# TODO(@LuckyMD): document relevant link at figshare.com/articles/*
URL = (
    "https://ndownloader.figshare.com/files/24566651?private_link=e3921450ec1bd0587870"
)


@utils.loader(data_url=URL, data_reference="https://doi.org/10.1126/science.aar4362")
def load_zebrafish(test=False):
    """Download zebrafish data from figshare."""
    with tempfile.TemporaryDirectory() as tempdir:
        filepath = os.path.join(tempdir, "zebrafish.h5ad")
        scprep.io.download.download_url(URL, filepath)
        adata = anndata.read_h5ad(filepath)
        utils.filter_genes_cells(adata)

    if test:
        adata = utils.subsample_even(adata, n_obs=500, even_obs="lab")
        utils.filter_genes_cells(adata)
        adata = adata[:, :100].copy()
        utils.filter_genes_cells(adata)
    return adata
