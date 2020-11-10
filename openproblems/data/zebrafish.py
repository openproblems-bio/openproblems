import os
import tempfile

import anndata
import scprep

from . import utils

URL = (
    "https://ndownloader.figshare.com/files/24566651?private_link=e3921450ec1bd0587870"
)


@utils.loader
def load_zebrafish(test=False):
    """Download zebrafish data from figshare."""
    with tempfile.TemporaryDirectory() as tempdir:
        filepath = os.path.join(tempdir, "zebrafish.h5ad")
        scprep.io.download.download_url(URL, filepath)
        adata = anndata.read_h5ad(filepath)
        adata.obs.index = adata.obs.index.astype(str)
        adata.obs.lab = adata.obs.lab.astype(str).astype("category")
        adata.obs.cell_type = adata.obs.cell_type.astype(str).astype("category")
        utils.filter_genes_cells(adata)

    if test:
        adata = utils.subsample_even(adata, n_obs=500, even_obs="lab")
        utils.filter_genes_cells(adata)
        adata = adata[:, :100]
        utils.filter_genes_cells(adata)
    return adata
