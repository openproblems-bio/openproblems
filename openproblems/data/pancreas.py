import os
import tempfile

import scprep
import scanpy as sc
import numpy as np

from .utils import loader, filter_genes_cells

URL = "https://ndownloader.figshare.com/files/24539828"


@loader
def load_pancreas(test=False):
    if test:
        # load full data first, cached if available
        adata = load_pancreas(test=False)
        import sys

        print(adata, file=sys.stderr)

        # Subsample pancreas data
        adata = adata[:, :500].copy()
        filter_genes_cells(adata)

        keep_celltypes = ["delta", "gamma"]
        keep_techs = ["inDrop4", "smarter", "celseq"]
        assert np.all(np.isin(keep_celltypes, adata.obs["celltype"])), adata.obs[
            "celltype"
        ].unique()
        assert np.all(np.isin(keep_techs, adata.obs["tech"])), adata.obs[
            "tech"
        ].unique()
        keep_tech_idx = adata.obs["tech"].isin(keep_techs)
        keep_celltype_idx = adata.obs["celltype"].isin(keep_celltypes)
        adata = adata[keep_tech_idx & keep_celltype_idx].copy()
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
    # Fix category dtypes
    adata.obs["tech"] = adata.obs["tech"].astype(str).astype("category")
    adata.obs["celltype"] = adata.obs["celltype"].astype(str).astype("category")

    # Remove processing
    adata.X = adata.layers["counts"]
    del adata.layers["counts"]
