from . import utils

import os
import scanpy as sc
import scprep
import tempfile

URL = "https://ndownloader.figshare.com/files/25717328"


@utils.loader
def load_immune(test=False):
    """Download immune human data from figshare."""
    if test:
        # load full data first, cached if available
        adata = load_immune(test=False)

        # Subsample pancreas data
        adata = adata[:, :500].copy()
        utils.filter_genes_cells(adata)

        # We are using the final annotation column from the original file
        # for the cell type annotation, the batch column contains the name
        # of the study the data originates from

        keep_final_annotations = adata.obs["final_annotation"].dtype.categories[[0, 3]]
        keep_batchs = adata.obs["batch"].dtype.categories[[0, -3, -2]]
        keep_batch_idx = adata.obs["batch"].isin(keep_batchs)
        keep_final_annotation_idx = adata.obs["final_annotation"].isin(
            keep_final_annotations
        )
        adata = adata[keep_batch_idx & keep_final_annotation_idx].copy()

        sc.pp.subsample(adata, n_obs=500)
        # Note: could also use 200-500 HVGs rather than 200 random genes

        # Ensure there are no cells or genes with 0 counts
        utils.filter_genes_cells(adata)

        return adata

    else:
        with tempfile.TemporaryDirectory() as tempdir:
            filepath = os.path.join(tempdir, "immune.h5ad")
            scprep.io.download.download_url(URL, filepath)
            adata = sc.read(filepath)

            # Note: anndata.X contains scran normalized data,
            # so we're storing it in layers['log_scran']
            adata.layers["log_scran"] = adata.X

            # Ensure there are no cells or genes with 0 counts
            utils.filter_genes_cells(adata)

        return adata
