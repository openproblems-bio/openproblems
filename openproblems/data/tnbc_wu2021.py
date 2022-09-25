from . import utils

import numpy as np
import os
import scanpy as sc
import scipy.sparse
import scprep
import tempfile

URL = "https://figshare.com/ndownloader/files/37593188"


@utils.loader(data_url=URL, data_reference="https://doi.org/10.1038/s41588-021-00911-1")
def load_tnbc_data(test=False):
    """Download TNBC data (Wu et al., 2021) from Figshare.

    Further information how we generated the assumed truth and the reference
    to the dataset is available at:
    https://figshare.com/articles/dataset/TNBC_Data_from_Wu_et_al_2021/20338536

    """
    if test:
        # load full data first, cached if available
        adata = load_tnbc_data(test=False)

        # Keep only relevant response ligands
        lr = list(set(adata.uns["ccc_target"].ligand))
        # add corresponding receptors
        lr += [
            "IFNLR1",
            "CD70",
            "NGFR",
            "CD53",
            "IL20RA",
            "IL20RB",
            "IL2RG",
            "SIGLEC6",
            "EGFR",
            "ITGB1",
            "ITGB3",
            "ACVRL1",
            "CXCR4",
            "SMAD3",
            "CAV1",
            "IL21R",
            "IL4R",
            "IL2RG",
            "IL13RA1",
            "IL13RA2",
        ]
        adata = adata[:, adata.var.index.isin(lr)].copy()

        # Keep only cells with nonzero relevant ligands/receptors
        utils.filter_genes_cells(adata)

        # Subsample data to 500 random cells
        sc.pp.subsample(adata, n_obs=500)

        # Add noise to ensure no genes are empty -- filtering genes could remove ligands
        empty_gene_idx = np.argwhere(adata.X.sum(axis=0).A.flatten() == 0).flatten()
        random_cell_idx = np.random.choice(adata.shape[0], len(empty_gene_idx))
        adata.X += scipy.sparse.coo_matrix(
            (
                np.ones(len(empty_gene_idx), dtype=adata.X.dtype),
                (random_cell_idx, empty_gene_idx),
            ),
            shape=adata.shape,
            dtype=adata.X.dtype,
        )
        adata.X = adata.X.tocsr()
        # update layers["counts"] with the new counts
        adata.layers["counts"] = adata.X
    else:
        with tempfile.TemporaryDirectory() as tempdir:
            filepath = os.path.join(tempdir, "brca_tnbc.h5ad")
            scprep.io.download.download_url(URL, filepath)
            adata = sc.read(filepath)

            # Ensure there are no cells or genes with 0 counts
            utils.filter_genes_cells(adata)

    # Remove unneeded cols
    adata.uns["ccc_target"] = adata.uns["ccc_target"][["ligand", "target", "response"]]

    return adata
