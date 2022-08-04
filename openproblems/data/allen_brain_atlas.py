from . import utils

import numpy as np
import os
import scanpy as sc
import scprep
import tempfile

URL = "https://figshare.com/ndownloader/files/36509385"


@utils.loader(data_url=URL, data_reference="https://doi.org/10.1038/nn.4216")
def load_mouse_brain_atlas(test=False):
    """Download Allen Brain (Taisc et al.,2016) data from Figshare."""
    if test:
        # load full data first, cached if available
        adata = load_mouse_brain_atlas(test=False)

        # Subsample data to 0.1 random cells
        sc.pp.subsample(adata, fraction=0.1)
        # Basic filtering
        utils.filter_genes_cells(adata)

        # Keep only 500 genes
        adata = adata[:, 30000:30500].copy()
        adata.uns["target_organism"] = adata.uns["target_organism"]

        # remove missing labels for tests
        labels = adata.obs["label"].cat.categories
        source_mask = ~adata.uns["ccc_target"]["target"].isin(labels)
        target_mask = ~adata.uns["ccc_target"]["target"].isin(labels)
        msk = np.logical_not(source_mask | target_mask)
        adata.uns["ccc_target"] = adata.uns["ccc_target"][msk]

        return adata

    else:
        with tempfile.TemporaryDirectory() as tempdir:
            filepath = os.path.join(tempdir, "mouse_brain_atlas.h5ad")
            scprep.io.download.download_url(URL, filepath)
            adata = sc.read(filepath)
            adata.uns["target_organism"] = adata.uns["target_organism"]

            # Ensure there are no cells or genes with 0 counts
            utils.filter_genes_cells(adata)

        return adata
