from .multimodal.scicar.cell_lines import rna_cells_url
from .multimodal.scicar.cell_lines import rna_genes_url
from .utils import loader

import anndata
import numpy as np
import pandas as pd
import scipy.sparse


@loader(data_url="https://openproblems.bio")
def load_sample_data(test=True):
    """Create a simple dataset to use for testing in multimodal applications."""
    assert test

    genes = pd.read_csv(rna_genes_url, low_memory=False, index_col=0).iloc[:500]
    cells = pd.read_csv(rna_cells_url, low_memory=False, index_col=0).iloc[:200]

    rna_data = scipy.sparse.csr_matrix(
        np.random.poisson(0.3, (cells.shape[0], genes.shape[0]))
    )

    adata = anndata.AnnData(rna_data, obs=cells, var=genes)
    adata.X = adata.X.astype(np.float64)
    return adata
