from .utils import loader

import anndata
import numpy as np


@loader
def load_template_data(test=False):
    if test:
        # Load the full dataset
        adata = load_template_data(test=False)
        # Subsample for speed
        adata = adata[:500, :200]
        return adata
    # Create or download the full dataset
    # TODO: update
    adata = anndata.AnnData(np.random.uniform(0, 1, size))
    return adata
