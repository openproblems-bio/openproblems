from . import utils

import os
import scprep
import tempfile

# sparsified from https://ndownloader.figshare.com/files/25555751
# TODO(@LuckyMD): change link to figshare.com/articles/*
URL = "https://ndownloader.figshare.com/files/36088649"


@utils.loader(data_url=URL, data_reference="nestorowa2016single")
def load_mouse_hspc_nestorowa2016(test=False):
    """Download Nestorowa data from Figshare."""
    import scanpy as sc

    if test:
        # load full data first, cached if available
        adata = load_mouse_hspc_nestorowa2016(test=False)

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
            filepath = os.path.join(tempdir, "human_blood_nestorowa2016.h5ad")
            scprep.io.download.download_url(URL, filepath)
            adata = sc.read(filepath)

            # Ensure there are no cells or genes with 0 counts
            utils.filter_genes_cells(adata)

        return adata
