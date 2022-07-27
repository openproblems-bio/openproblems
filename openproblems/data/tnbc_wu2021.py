from . import utils

import os
import scanpy as sc
import scprep
import tempfile

URL = "https://figshare.com/ndownloader/files/36352914"


@utils.loader(data_url=URL, data_reference="https://doi.org/10.1038/s41588-021-00911-1")
def load_tnbc_data(test=False):
    """Download TNBC data (Wu et al., 2021) from Figshare."""
    if test:
        # load full data first, cached if available
        adata = load_tnbc_data(test=False)

        # Keep only relevant response ligands
        lr = list(set(adata.uns["bench"].ligand))
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

        # Remove empty
        utils.filter_genes_cells(adata)

        return adata

    else:
        with tempfile.TemporaryDirectory() as tempdir:
            filepath = os.path.join(tempdir, "brca_tnbc.h5ad")
            scprep.io.download.download_url(URL, filepath)
            adata = sc.read(filepath)

            # Ensure there are no cells or genes with 0 counts
            utils.filter_genes_cells(adata)

        return adata
