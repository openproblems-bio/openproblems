from .....data.immune_cells import load_immune
from .....tools.decorators import dataset

import scanpy as sc


@dataset(dataset_name="Immune (by batch)", image="openproblems")
def immune_batch(test=False):
    adata = load_immune(test)
    from_cache = adata.__from_cache__
    adata.obs["labels"] = adata.obs["final_annotation"]

    sc.pp.filter_genes(adata, min_counts=1)
    sc.pp.filter_genes(adata, min_cells=1)

    sc.tl.pca(
        adata,
        svd_solver="arpack",
        return_info=True,
    )
    adata.obsm["X_uni"] = adata.obsm["X_pca"]

    sc.pp.neighbors(adata, use_rep="X_uni", key_added="uni")

    adata.__from_cache__ = from_cache
    return adata
