import os
import tempfile

import scprep
import scanpy as sc

from .utils import loader, filter_genes_cells

URL = "https://ndownloader.figshare.com/files/25555751"


@loader
def load_nestorowa_2016(test=False):
    if test:
        # load full data first, cached if available
        adata = load_nestorowa_2016(test=False)

        # Subsample pancreas data
        adata = adata[:, :500].copy()
        filter_genes_cells(adata)

        sc.pp.subsample(adata, n_obs=500)
        # Note: could also use 200-500 HVGs rather than 200 random genes

        # Ensure there are no cells or genes with 0 counts
        filter_genes_cells(adata)

        return adata

    else:
        with tempfile.TemporaryDirectory() as tempdir:
            filepath = os.path.join(tempdir, "nestorowa_2016.h5ad")
            scprep.io.download.download_url(URL, filepath)
            adata = sc.read(filepath)

            # Ensure there are no cells or genes with 0 counts
            filter_genes_cells(adata)

        return adata

