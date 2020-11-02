import os
import tempfile

import scprep
import scanpy as sc

from .utils import loader, filter_genes_cells

URL = "https://ndownloader.figshare.com/files/24539828"


@loader
def load_pancreas(test=False):
    if test:
        # load full data first, cached if available
        adata = load_pancreas(test=False)

        # Subsample pancreas data
        cts = ["delta", "gamma"]
        batches = ["inDrop4", "smarter", "celseq"]
        adata = adata[
            (adata.obs["batch"].isin(batches)) & (adata.obs["labels"].isin(cts)), :500
        ].copy()
        # Note: could also use 200-500 HVGs rather than 200 random genes

        # Ensure there are no cells or genes with 0 counts
        filter_genes_cells(adata)

        return adata

    else:
        with tempfile.TemporaryDirectory() as tempdir:
            filepath = os.path.join(tempdir, "pancreas.h5ad")
            scprep.io.download.download_url(URL, filepath)
            adata = sc.read(filepath)

            # Correct covariate naming and remove preprocessing
            prep_pancreas(adata)

            # Ensure there are no cells or genes with 0 counts
            filter_genes_cells(adata)

        return adata


def prep_pancreas(adata):
    # Rename categories
    adata.obs["batch"] = adata.obs["tech"].copy()
    adata.obs["labels"] = adata.obs["celltype"].copy()

    # Drop excess covariates
    adata.obs = adata.obs.drop(columns=["tech", "size_factors", "celltype"])

    # Remove processing
    adata.X = adata.layers["counts"]
    del adata.layers["counts"]
