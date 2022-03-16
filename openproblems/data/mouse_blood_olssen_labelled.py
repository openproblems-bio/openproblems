from . import utils

import os
import scanpy as sc
import scprep
import tempfile

URL = "https://ndownloader.figshare.com/files/27346712"


@utils.loader
def load_olsson_2016_mouse_blood(test=False):
    """Download Olsson, 2016_mouse_blood, Nature, 2016 data from Figshare."""
    if test:
        # load full data first, cached if available
        adata = load_olsson_2016_mouse_blood(test=False)

        # Subsample data
        adata = adata[:, :500].copy()
        utils.filter_genes_cells(adata)

        sc.pp.subsample(adata, n_obs=500)
        # Note: could also use 200-500 HVGs rather than 200 random genes

        # Ensure there are no cells or genes with 0 counts
        utils.filter_genes_cells(adata)

        return adata

    else:
        with tempfile.TemporaryDirectory() as tempdir:
            filepath = os.path.join(tempdir, "olsson_2016_mouse_blood.h5ad")
            scprep.io.download.download_url(URL, filepath)
            adata = sc.read(filepath)

            # Ensure there are no cells or genes with 0 counts
            utils.filter_genes_cells(adata)

        return adata
