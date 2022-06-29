from . import utils

import os
import scanpy as sc
import scprep
import tempfile

URL = "https://ndownloader.figshare.com/files/36088646"


@utils.loader(data_url=URL)
def load_olsson_2016_mouse_blood(test=False):
    """Download Olsson, 2016_mouse_blood, Nature, 2016 data from Figshare."""
    if test:
        # load full data first, cached if available
        adata = load_olsson_2016_mouse_blood(test=False)

        # Select 700 genes expressed on most cells
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        top = list(adata.var.nlargest(700, "mean_counts").index.values)
        adata = adata[:, top].copy()
        utils.filter_genes_cells(adata)

        # Subsample to 300 cells
        sc.pp.subsample(adata, n_obs=300)
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
