import pandas as pd
import scprep
import tempfile
import os

from ..utils import filter_joint_data_empty_cells, create_joint_adata


def load_scicar(
    rna_url,
    rna_cells_url,
    rna_genes_url,
    atac_url,
    atac_cells_url,
    atac_genes_url,
):
    rna_cells = pd.read_csv(rna_cells_url, low_memory=False)["sample"]
    rna_genes = pd.read_csv(rna_genes_url, low_memory=False)["gene_id"]
    atac_cells = pd.read_csv(atac_cells_url, low_memory=False)["sample"]
    atac_genes = pd.read_csv(atac_genes_url, low_memory=False, index_col=0)["peak"]

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
        X_index=rna_cells,
        X_columns=rna_genes,
        Y_index=atac_cells,
        Y_columns=atac_genes,
    )
    adata = filter_joint_data_empty_cells(adata)
    return adata
