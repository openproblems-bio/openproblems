from .....data.pancreas import load_pancreas
from .....tools.decorators import dataset

import scanpy as sc


@dataset(
    dataset_name="Pancreas (by batch)",
    data_url=load_pancreas.metadata["data_url"],
    dataset_summary="TODO",
    image="openproblems",
)
def pancreas_batch(test=False):
    adata = load_pancreas(test)
    adata.obs["labels"] = adata.obs["celltype"]
    adata.obs["batch"] = adata.obs["tech"]

    sc.pp.filter_genes(adata, min_counts=1)
    sc.pp.filter_genes(adata, min_cells=1)

    sc.tl.pca(
        adata,
        svd_solver="arpack",
        return_info=True,
    )
    adata.obsm["X_uni"] = adata.obsm["X_pca"]

    sc.pp.neighbors(adata, use_rep="X_uni", key_added="uni")

    return adata
