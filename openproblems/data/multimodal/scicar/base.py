import pandas as pd
import scprep
import tempfile
import os
import numpy as np
import anndata

from ..utils import filter_joint_data_empty_cells, create_joint_adata


def load_scicar(
    rna_url,
    rna_cells_url,
    rna_genes_url,
    atac_url,
    atac_cells_url,
    atac_genes_url,
    with_peak=False,
):
    if not with_peak:
        rna_cells = pd.read_csv(rna_cells_url, low_memory=False)["sample"]
        rna_genes = pd.read_csv(rna_genes_url, low_memory=False)["gene_id"]
        atac_cells = pd.read_csv(atac_cells_url, low_memory=False)["sample"]
        atac_genes = pd.read_csv(atac_genes_url, low_memory=False, index_col=0)["peak"]
    else:
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

    print(atac_genes)

    if not with_peak:
        adata = create_joint_adata(
            rna_data,
            atac_data,
            X_index=rna_cells,
            X_columns=rna_genes,
            Y_index=atac_cells,
            Y_columns=atac_genes,
        )
        adata = filter_joint_data_empty_cells(adata)
    else:
        rna_data, rna_cells_index = scprep.filter.filter_empty_cells(
            rna_data, rna_cells.index
        )
        atac_data, atac_cells_index = scprep.filter.filter_empty_cells(
            atac_data, atac_cells.index
        )

        common_cells = np.intersect1d(rna_cells_index, atac_cells_index)

        rna_subset = np.isin(rna_cells_index, common_cells)
        rna_cells, rna_data = rna_cells.loc[rna_subset], rna_data[rna_subset]
        rna_order = np.argsort(rna_cells.index)
        rna_cells, rna_data = rna_cells.iloc[rna_order], rna_data[rna_order]

        atac_subset = np.isin(atac_cells_index, common_cells)
        atac_cells, atac_data = atac_cells.loc[atac_subset], atac_data[atac_subset]
        atac_order = np.argsort(atac_cells.index)
        atac_cells, atac_data = atac_cells.iloc[atac_order], atac_data[atac_order]

        adata = anndata.AnnData(
            rna_data,
            obs=rna_cells,
            var=rna_genes,
        )
        adata.obsm["mode2"] = atac_data
        adata.uns["mode2_obs"] = atac_cells.to_numpy()
        adata.uns["mode2_var"] = atac_genes.iloc[:, 1:].to_numpy()
    return adata
