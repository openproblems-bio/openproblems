from ..utils import create_joint_adata
from ..utils import filter_joint_data_empty_cells

import numpy as np
import os
import pandas as pd
import scprep
import tempfile


def load_snare(
    rna_url,
    rna_cells_url,
    rna_genes_url,
    atac_url,
    atac_cells_url,
    atac_genes_url,
):
    """Load sci-CAR data from GEO."""
    rna_genes = pd.read_csv(
        rna_genes_url, low_memory=False, index_col=0, compression="gzip", header=None
    )
    atac_genes = pd.read_csv(
        atac_genes_url, low_memory=False, index_col=0, compression="gzip", header=None
    )
    rna_cells = pd.read_csv(
        rna_cells_url, low_memory=False, index_col=0, compression="gzip", header=None
    )
    atac_cells = pd.read_csv(
        atac_cells_url, low_memory=False, index_col=0, compression="gzip", header=None
    )

    with tempfile.TemporaryDirectory() as tempdir:
        rna_file = os.path.join(tempdir, "rna.mtx.gz")
        scprep.io.download.download_url(rna_url, rna_file)
        rna_data = scprep.io.load_mtx(rna_file, cell_axis="col").tocsr()
        atac_file = os.path.join(tempdir, "atac.mtx.gz")
        scprep.io.download.download_url(atac_url, atac_file)
        atac_data = scprep.io.load_mtx(atac_file, cell_axis="col").tocsr()

    adata = create_joint_adata(
        rna_data,
        atac_data,
        X_index=rna_cells.index,
        X_columns=rna_genes.index,
        Y_index=atac_cells.index,
        Y_columns=atac_genes.index,
    )

    # merge obs and var
    adata.obs = rna_cells.loc[adata.obs.index]
    rna_genes.loc[:, "gene_short_name"] = rna_genes.index
    rna_genes.index.name = "genes"
    adata.var = rna_genes

    adata.uns["mode2_var"] = np.array(
        list(
            map(
                lambda x: [x.split(":")[0]] + x.split(":")[1].split("-"),
                list(adata.uns["mode2_var"]),
            )
        )
    )
    adata.uns["mode2_var_chr"] = adata.uns["mode2_var"][:, 0]
    adata.uns["mode2_var_start"] = adata.uns["mode2_var"][:, 1]
    adata.uns["mode2_var_end"] = adata.uns["mode2_var"][:, 2]
    adata = filter_joint_data_empty_cells(adata)
    return adata
