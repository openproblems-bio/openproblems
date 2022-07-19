from . import utils

import os
import scanpy as sc
import scprep
import tempfile

URL = "https://ndownloader.figshare.com/files/25717328"


@utils.loader(data_url=URL, data_reference="TODO")
def load_immune(test=False):
    """Download immune human data from figshare."""
    if test:
        # load full data first, cached if available
        adata = load_immune(test=False)

        # Subsample immune data to two batches with 250 cells each
        adata = adata[:, :500].copy()
        batch1 = adata[adata.obs.batch == "Oetjen_A"][:250]
        batch2 = adata[adata.obs.batch == "Freytag"][:250]
        adata = batch1.concatenate(batch2)
        # Note: could also use 200-500 HVGs rather than 200 random genes

        # Ensure there are no cells or genes with 0 counts
        utils.filter_genes_cells(adata)

        return adata

    else:
        with tempfile.TemporaryDirectory() as tempdir:
            filepath = os.path.join(tempdir, "immune.h5ad")
            scprep.io.download.download_url(URL, filepath)
            adata = sc.read(filepath)

            # Note: anndata.X contains scran log-normalized data,
            # so we're storing it in layers['log_scran']
            adata.layers["log_scran"] = adata.X
            adata.X = adata.layers["counts"]
            del adata.layers["counts"]

            # Ensure there are no cells or genes with 0 counts
            utils.filter_genes_cells(adata)

        return adata
