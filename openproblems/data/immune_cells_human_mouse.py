from . import utils

import os
import scprep
import tempfile

# sparsified from https://figshare.com/articles/dataset/Benchmarking_atlas-level_data_integration_in_single-cell_genomics_-_integration_task_datasets_Immune_and_pancreas_/12420968/2 # noqa: E501
URL = "https://figshare.com/ndownloader/files/25717166"


@utils.loader(data_url=URL, data_reference="luecken2022benchmarking")
def load_immune_hm(test=False):
    """Download immune human data from figshare."""
    import scanpy as sc

    if test:
        # load full data first, cached if available
        adata = load_immune_hm(test=False)

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

            # NOTE: adata.X contains log-normalized data, so we're moving it
            adata.layers["log_normalized"] = adata.X
            adata.X = adata.layers["counts"]

            # Ensure there are no cells or genes with 0 counts
            utils.filter_genes_cells(adata)

        return adata
