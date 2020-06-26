import numpy as np
import anndata


def dummy(test=False):
    adata = anndata.AnnData(np.random.uniform(0, 1, (100, 10)))
    adata.obs["labels"] = np.random.choice(2, 100, replace=True)
    adata.obs["is_train"] = np.random.choice([True, False], 100, replace=True)
    return adata
