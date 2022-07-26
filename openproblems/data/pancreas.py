from . import utils

import anndata as ad
import numpy as np
import os
import scanpy as sc
import scprep
import tempfile

URL = "https://ndownloader.figshare.com/files/24539828"


@utils.loader(data_url=URL, data_reference="https://doi.org/10.1038/s41592-021-01336-8")
def load_pancreas(test=False, integer_only=False):
    """Download pancreas data from figshare."""
    if test:
        # load full data first, cached if available
        adata = load_pancreas(test=False, integer_only=integer_only)

        keep_celltypes = adata.obs["celltype"].dtype.categories[[0, 3]]
        keep_techs = adata.obs["tech"].dtype.categories[[0, -3, -2]]
        keep_tech_idx = adata.obs["tech"].isin(keep_techs)
        keep_celltype_idx = adata.obs["celltype"].isin(keep_celltypes)
        adata = adata[keep_tech_idx & keep_celltype_idx].copy()

        # Subsample pancreas data
        adata = adata[:, :500].copy()
        # Note: could also use 200-500 HVGs rather than 200 random genes
        utils.filter_genes_cells(adata)

        # select 250 cells from each celltype
        keep_cell_idx = np.concatenate(
            [
                np.random.choice(
                    np.argwhere(adata.obs["celltype"].to_numpy() == ct).flatten(),
                    250,
                    replace=False,
                )
                for ct in keep_celltypes
            ]
        )
        adata = adata[adata.obs[keep_cell_idx]]

        # Ensure there are no cells or genes with 0 counts
        utils.filter_genes_cells(adata)

        return adata

    with tempfile.TemporaryDirectory() as tempdir:
        filepath = os.path.join(tempdir, "pancreas.h5ad")
        scprep.io.download.download_url(URL, filepath)
        adata = sc.read(filepath)

    # NOTE: X contains counts that are normalized with scran
    adata.layers["log_scran"] = adata.X
    adata.X = adata.layers["counts"]
    del adata.layers["counts"]

    if integer_only:
        adata = _get_pancreas_integer(adata)

    # Ensure there are no cells or genes with 0 counts
    utils.filter_genes_cells(adata)

    return adata


def _get_pancreas_integer(adata: ad.AnnData):
    """Transform counts to integer.

    For some platforms the pancreas data set only have processed counts.
    Here we grab those with integer counts.
    See https://github.com/theislab/scib-reproducibility/tree/main/notebooks/data_preprocessing/pancreas # noqa: E501
    """
    is_int = ["smartseq2"]
    is_int += ["inDrop{}".format(x) for x in range(1, 5)]

    keep = np.zeros(len(adata)).astype(bool)

    for tech in is_int:
        idx = adata.obs.tech.values == tech
        keep = keep | idx

    adata = adata[keep, :].copy()

    return adata
