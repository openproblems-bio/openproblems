import os
import tempfile

import scprep
import scanpy as sc

from .utils import loader, filter_genes_cells

URL = "https://ndownloader.figshare.com/files/24539828"


@loader
def load_pancreas(test=False):
    """Download pancreas data from figshare."""
    if test:
        # load full data first, cached if available
        adata = load_pancreas(test=False)

        # Subsample pancreas data
        adata = adata[:, :500].copy()
        filter_genes_cells(adata)

        keep_celltypes = adata.obs["celltype"].dtype.categories[[0, 3]]
        keep_techs = adata.obs["tech"].dtype.categories[[0, -3, -2]]
        keep_tech_idx = adata.obs["tech"].isin(keep_techs)
        keep_celltype_idx = adata.obs["celltype"].isin(keep_celltypes)
        adata = adata[keep_tech_idx & keep_celltype_idx].copy()

        sc.pp.subsample(adata, n_obs=500)
        # Note: could also use 200-500 HVGs rather than 200 random genes

        # Ensure there are no cells or genes with 0 counts
        filter_genes_cells(adata)

        return adata

    else:
        with tempfile.TemporaryDirectory() as tempdir:
            filepath = os.path.join(tempdir, "pancreas.h5ad")
            scprep.io.download.download_url(URL, filepath)
            adata = sc.read(filepath)

            # Remove preprocessing
            adata.X = adata.layers["counts"]
            del adata.layers["counts"]

            # Ensure there are no cells or genes with 0 counts
            filter_genes_cells(adata)

        return adata
