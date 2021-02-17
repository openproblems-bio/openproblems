from ..utils import create_joint_adata
from ..utils import filter_joint_data_empty_cells
from scipy import sparse

import numpy as np
import os


def load_share(
    rna_url,
    atac_url,
    atac_cells_url,
    atac_genes_url,
):
    """Load SHARE-seq data from GEO."""
    import datatable as dt
    import pandas as pd
    import scprep
    import tempfile

    with tempfile.TemporaryDirectory() as tempdir:
        atac_gene_file = os.path.join(tempdir, "atac_gene.csv.gz")
        scprep.io.download.download_url(atac_genes_url, atac_gene_file)
        atac_genes = pd.read_table(atac_gene_file, header=None)

        atac_cell_file = os.path.join(tempdir, "atac_cell.csv.gz")
        scprep.io.download.download_url(atac_cells_url, atac_cell_file)
        atac_cells = pd.read_table(atac_cell_file, header=None)

        rna_file = os.path.join(tempdir, "rna.csv.gz")
        scprep.io.download.download_url(rna_url, rna_file)
        rna_data = dt.fread(rna_file)
        rna_cells = np.array(list(map(lambda x: x.replace(",", "."), rna_data.keys())))[
            1:
        ]
        rna_genes = np.array(rna_data[:, "gene"].to_list()[0])
        rna_data = rna_data[:, 1:].to_numpy().transpose()
        rna_data = sparse.csr_matrix(rna_data)

        atac_file = os.path.join(tempdir, "atac.mtx.gz")
        scprep.io.download.download_url(atac_url, atac_file)
        atac_data = scprep.io.load_mtx(atac_file, cell_axis="col").tocsr()
        atac_cells = np.array(
            list(map(lambda x: x.replace(",", "."), atac_cells.iloc[:, 0].values))
        )

        print(rna_data.shape)
        print(atac_data.shape)

    adata = create_joint_adata(
        rna_data,
        atac_data,
        X_index=rna_cells,
        X_columns=rna_genes,
        Y_index=atac_cells,
        Y_columns=atac_genes.iloc[:, 0:3].values,
    )

    adata.uns["mode2_var_chr"] = adata.uns["mode2_var"][:, 0]
    adata.uns["mode2_var_start"] = adata.uns["mode2_var"][:, 1]
    adata.uns["mode2_var_end"] = adata.uns["mode2_var"][:, 2]

    adata.var["gene_short_name"] = adata.var.index
    adata = filter_joint_data_empty_cells(adata)
    return adata
