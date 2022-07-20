from . import utils

import os
import scanpy as sc
import scprep
import tempfile

URL = "https://figshare.com/ndownloader/files/36352287"


@utils.loader(data_url=URL)
def load_mouse_brain_atlas(test=False):
    """Download Allen Brain (Taisc et al.,2016) data from Figshare."""
    if test:
        # load full data first, cached if available
        adata = load_mouse_brain_atlas(test=False)

        # Subsample data to 500 random cells
        sc.pp.subsample(adata, n_obs=500)
        # Basic filtering
        utils.filter_genes_cells(adata)

        # Keep only 500 genes
        adata = adata[:, 30000:30500].copy()
        adata.uns["target_organism"] = adata.uns["target_organism"][0]

        return adata

    else:
        with tempfile.TemporaryDirectory() as tempdir:
            filepath = os.path.join(tempdir, "mouse_brain_atlas.h5ad")
            scprep.io.download.download_url(URL, filepath)
            adata = sc.read(filepath)
            adata.uns["target_organism"] = adata.uns["target_organism"][0]

            # Ensure there are no cells or genes with 0 counts
            utils.filter_genes_cells(adata)

        return adata
