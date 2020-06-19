import numpy as np
import anndata


def dummy():
    adata = anndata.AnnData(np.random.uniform(0, 1, (100, 10)))
    adata.obs["labels"] = np.random.choice(2, 100, replace=True)
    return adata
