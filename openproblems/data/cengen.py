from . import utils

import os
import scanpy as sc
import scprep
import tempfile

URL = (
    "https://github.com/Munfred/wormcells-data/"
    "releases/download/taylor2020/taylor2020.h5ad"
)


@utils.loader
def load_cengen(test=False):
    """Download CeNGEN data from GitHub.

    To learn about how the data was generated visit www.cengen.org
    To learn about WormBase curation efforts for C. elegans single cell
    data visit https://wormbase.github.io/single-cell/
    """
    with tempfile.TemporaryDirectory() as tempdir:
        filepath = os.path.join(tempdir, "cengen.h5ad")
        scprep.io.download.download_url(URL, filepath)
        adata = sc.read_h5ad(filepath)
        utils.filter_genes_cells(adata)

    if test:
        batches = adata.obs.experiment_code.cat.categories[3:6]
        adata = adata[adata.obs.experiment_code.isin(batches), :]
        sc.pp.subsample(adata, n_obs=500)
        adata = adata[:, :1000]
        utils.filter_genes_cells(adata)

    return adata
