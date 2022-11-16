from . import utils

import os
import scprep
import tempfile

# TODO(@LuckyMD): document relevant link at figshare.com/articles/*
URL = (
    "https://github.com/Munfred/wormcells-data/"
    "releases/download/taylor2020/taylor2020.h5ad"
)


@utils.loader(
    data_url=URL, data_reference="https://doi.org/10.1016/j.neuron.2018.07.042"
)
def load_cengen(test=False):
    """Download CeNGEN data from GitHub.

    To learn about how the data was generated visit www.cengen.org
    To learn about WormBase curation efforts for C. elegans single cell
    data visit https://wormbase.github.io/single-cell/
    """
    import scanpy as sc

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
