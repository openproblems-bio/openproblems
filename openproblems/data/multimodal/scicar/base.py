from ..utils import create_joint_adata
from ..utils import filter_joint_data_empty_cells

import os
import pandas as pd
import scprep
import tempfile

DATA_REFERENCE = "https://doi.org/10.1126/science.aau0730"


def load_scicar(
    rna_url,
    rna_cells_url,
    rna_genes_url,
    atac_url,
    atac_cells_url,
    atac_genes_url,
):
    """Load sci-CAR data from GEO."""
    rna_genes = pd.read_csv(rna_genes_url, low_memory=False, index_col=0)
    atac_genes = pd.read_csv(atac_genes_url, low_memory=False, index_col=1)
    rna_cells = pd.read_csv(rna_cells_url, low_memory=False, index_col=0)
    atac_cells = pd.read_csv(atac_cells_url, low_memory=False, index_col=0)

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
    adata.var = rna_genes
    for key in atac_cells.columns:
        adata.obs[key] = atac_cells[key]
    adata.uns["mode2_varnames"] = []
    for key in atac_genes.columns:
        varname = "mode2_var_{}".format(key)
        adata.uns[varname] = atac_genes[key].values
        adata.uns["mode2_varnames"].append(varname)

    adata = filter_joint_data_empty_cells(adata)
    return adata
