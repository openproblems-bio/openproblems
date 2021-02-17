from ....data.pancreas import load_pancreas
from ....tools.decorators import dataset
from scIB.preprocessing import normalize, hvg_batch

import numpy as np
import scanpy as sc


@dataset("Pancreas (by batch)")
def pancreas_batch(test=False):
    adata = load_pancreas(test=test)
    adata.obs["labels"] = adata.obs["celltype"]

    adata.obs["batch"] = adata.obs["tech"]

    normalize(adata)

    hvg_list = hvg_batch(adata, batch_key="batch")
    adata.var["highly_variable"] = np.in1d(adata.var_names, hvg_list)
    sc.tl.pca(
        adata,
        use_highly_variable="highly_variable",
        svd_solver="arpack",
        return_info=True,
    )
    adata.obsm["X_uni"] = adata.obsm["X_pca"]

    sc.pp.neighbors(adata, use_rep="X_uni", key_added="uni")

    return adata
