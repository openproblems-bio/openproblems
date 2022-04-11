from . import utils

import os
import scanpy as sc
import scprep
import tempfile

URL = "https://ndownloader.figshare.com/files/25555739"


@utils.loader
def load_tenx_5k_pbmc(test=False):
    """Download 5k PBMCs from 10x Genomics."""
    if test:
        # load full data first, cached if available
        adata = load_tenx_5k_pbmc(test=False)

        # Subsample pancreas data
        adata = adata[:, :500].copy()
        utils.filter_genes_cells(adata)

        sc.pp.subsample(adata, n_obs=500)
        # Note: could also use 200-500 HVGs rather than 200 random genes

        # Ensure there are no cells or genes with 0 counts
        utils.filter_genes_cells(adata)

        return adata

    else:
        with tempfile.TemporaryDirectory() as tempdir:
            filepath = os.path.join(tempdir, "10x_5k_pbmc.h5ad")
            scprep.io.download.download_url(URL, filepath)
            adata = sc.read(filepath)

            adata.var_names_make_unique()

            # Ensure there are no cells or genes with 0 counts
            utils.filter_genes_cells(adata)

        return adata
