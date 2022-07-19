from .....data.immune_cells import load_immune
from .....tools.decorators import dataset

import scanpy as sc


@dataset(
    dataset_name="Immune (by batch)",
    data_url=load_immune.metadata["data_url"],
    dataset_summary="Human immune cells from peripheral blood and bone marrow "
    "taken from 5 datasets comprising 10 batches across technologies (10X, "
    "Smart-seq2). Data was retrieved from doi:10.1038/s41592-021-01336-8.",
    image="openproblems",
)
def immune_batch(test=False):
    adata = load_immune(test)
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
    return adata
