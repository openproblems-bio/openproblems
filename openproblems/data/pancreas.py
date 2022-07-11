from . import utils

import anndata as ad
import numpy as np
import os
import scanpy as sc
import scprep
import tempfile

URL = "https://ndownloader.figshare.com/files/24539828"


@utils.loader(data_url=URL)
def load_pancreas(test=False):
    """Download pancreas data from figshare."""
    if test:
        # load full data first, cached if available
        adata = load_pancreas(test=False)

        # Subsample pancreas data
        adata = adata[:, :500].copy()
        utils.filter_genes_cells(adata)

        keep_celltypes = adata.obs["celltype"].dtype.categories[[0, 3]]
        keep_techs = adata.obs["tech"].dtype.categories[[0, -3, -2]]
        keep_tech_idx = adata.obs["tech"].isin(keep_techs)
        keep_celltype_idx = adata.obs["celltype"].isin(keep_celltypes)
        adata = adata[keep_tech_idx & keep_celltype_idx].copy()

        sc.pp.subsample(adata, n_obs=500)
        # Note: could also use 200-500 HVGs rather than 200 random genes

        # Ensure there are no cells or genes with 0 counts
        utils.filter_genes_cells(adata)

        return adata

    else:
        with tempfile.TemporaryDirectory() as tempdir:
            filepath = os.path.join(tempdir, "pancreas.h5ad")
            scprep.io.download.download_url(URL, filepath)
            adata = sc.read(filepath)

            # NOTE: X contains counts that are normalized with scran
            adata.layers["log_scran"] = adata.X
            adata.X = adata.layers["counts"]
            del adata.layers["counts"]
            # Ensure there are no cells or genes with 0 counts
            utils.filter_genes_cells(adata)

        return adata


def get_pancreas_integer(adata: ad.AnnData):
    """Transform counts to integer.

    For some platforms the pancreas data set only have processed counts.
    Here we grab those with integer counts.
    See: https://github.com/theislab/scib-reproducibility/tree/main/notebooks/data_preprocessing/pancreas. # noqa: E501
    """
    is_int = ["smartseq2"]
    is_int += ["inDrop{}".format(x) for x in range(1, 5)]

    keep = np.zeros(len(adata)).astype(bool)

    for tech in is_int:
        idx = adata.obs.tech.values == tech
        keep = keep | idx

    adata = adata[keep, :].copy()

    return adata
