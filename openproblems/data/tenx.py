from . import utils

import os
import scanpy as sc
import scprep
import tempfile

PBMC_1K_URL = "https://ndownloader.figshare.com/files/24974582"
PBMC_5K_URL = "https://ndownloader.figshare.com/files/25555739"
REFERENCE_URL = "https://www.10xgenomics.com/resources/datasets"


@utils.loader(data_url=PBMC_1K_URL, data_reference=REFERENCE_URL)
def load_tenx_1k_pbmc(test=False):
    """Download PBMC data from Figshare."""
    if test:
        adata = load_tenx_1k_pbmc(test=False)
        sc.pp.subsample(adata, n_obs=100)
        adata = adata[:, :1000]
        utils.filter_genes_cells(adata)
    else:
        with tempfile.TemporaryDirectory() as tempdir:
            filepath = os.path.join(tempdir, "pbmc.h5ad")
            scprep.io.download.download_url(PBMC_1K_URL, filepath)
            adata = sc.read_h5ad(filepath)
            utils.filter_genes_cells(adata)
    return adata


@utils.loader(data_url=PBMC_5K_URL, data_reference=REFERENCE_URL)
def load_tenx_5k_pbmc(test=False):
    """Download 5k PBMCs from 10x Genomics."""
    if test:
        # load full data first, cached if available
        adata = load_tenx_5k_pbmc(test=False)

        # Subsample pancreas data
        adata = adata[:, :500].copy()
        utils.filter_genes_cells(adata)

        sc.pp.subsample(adata, n_obs=500)
        # Note: could also use 200-500 HVGs rather than 200 random genes

        # Ensure there are no cells or genes with 0 counts
        utils.filter_genes_cells(adata)

        return adata

    else:
        with tempfile.TemporaryDirectory() as tempdir:
            filepath = os.path.join(tempdir, "10x_5k_pbmc.h5ad")
            scprep.io.download.download_url(PBMC_5K_URL, filepath)
            adata = sc.read(filepath)

            adata.var_names_make_unique()

            # Ensure there are no cells or genes with 0 counts
            utils.filter_genes_cells(adata)

        return adata
