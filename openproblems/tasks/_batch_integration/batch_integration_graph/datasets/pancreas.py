from .....data.pancreas import load_pancreas
from .....tools.decorators import dataset
from .....tools.normalize import log_scran_pooling
#from scIB.preprocessing import normalize, hvg_batch

import numpy as np
import scanpy as sc


@dataset("Pancreas (by batch)")
def pancreas_batch(test=False):
    adata = load_pancreas(test=test)
    adata.obs["labels"] = adata.obs["celltype"]

    adata.obs["batch"] = adata.obs["tech"]

    log_scran_pooling(adata)
    adata.layers['normalized'] = adata.X

    sc.tl.pca(
        adata,
        svd_solver="arpack",
        return_info=True,
    )
    adata.obsm["X_uni"] = adata.obsm["X_pca"]

    sc.pp.neighbors(adata, use_rep="X_uni", key_added="uni")

    if test:
        sc.pp.subsample(adata, n_obs=200)
        return adata[:, :500]
    return adata
